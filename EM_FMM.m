function [VBN_EM, SumError] = EM_FMM(pixels,Init,upNum,distribution,Flag,address_fix)
% Flag=1:显示参数迭代更新图figure(3)，否则不显示
% upNum:迭代上限数值
% Init：迭代之前VBN的初始参数集，格式为
% [Mu Var W
%  ↓ ↓  ↓]
% pixels 就是输入的图的像素向量
% distribution是表示分布的字符数组，比如，如果是三高斯，distribution='ggg'，如果是两指数和一高斯，distribution='eeg'
% address_fix 指固定不更新的分布的索引，比如不更新第4,5个分布，那么fix = [4,5]

LengthRegion = numel(pixels);

% 根据distribution字符串计算分布的位置顺序
address_gaussian = find(distribution == 'g');
address_rayleigh = find(distribution == 'r');
address_exponent = find(distribution == 'e');
% 固定分布默认没有
if nargin == 5
    address_fix = [];
end

K = size(Init,1);
N = length(pixels(:));
% 如果有e分布，一定要从1开始，从1开始，从1开始！
[N_pixels,X_pixels] = hist(pixels(:),min(pixels):max(pixels));
L = length(min(pixels):max(pixels));
[N_pixels,X_pixels] = hist(pixels(:),min(pixels):max(pixels));
L = length(min(pixels):max(pixels));

FC = M_Extand(N_pixels,K,L);
XI = M_Extand(X_pixels,K,L);

Mu = zeros(K,upNum+1);
Var = zeros(K,upNum+1);
W = zeros(K,upNum+1);
Error = zeros(1,upNum);
dL = length(Error);

Mu(:,1) = Init(:,1);
Var(:,1) = Init(:,2).^2;
W(:,1) = Init(:,3);

for i = 1:upNum
    [plx,pxl] = pLX(min(pixels):max(pixels),K,Mu(:,i),Var(:,i),W(:,i),address_gaussian,address_rayleigh,address_exponent);
    W(:,i+1) = (1/N)*sum(FC.*plx,2);
%     Mu(1:K-1,i+1) = sum(FC(1:K-1,:).*plx(1:K-1,:),2)./(sum(XI(1:K-1,:).*FC(1:K-1,:).*plx(1:K-1,:),2));% (1:K-1,:)对应指数分布密度函数
%     Mu(K,i+1) = sum(XI(K,:).*FC(K,:).*plx(K,:),2)./sum(FC(K,:).*plx(K,:),2);% 区间(K,:)对应的高斯均值
%     %Mu(K,i+1) = Mu(K);% 区间K对应的高斯均值
%     MU = M_Extand(Mu(:,i+1),K,L);
%     Var(K,i+1) = sum((XI(K,:)-MU(K,:)).^2.*FC(K,:).*plx(K,:),2)./sum(FC(K,:).*plx(K,:),2);
    % 对应gaussian分布均值
    Mu(address_gaussian, i+1) = sum(XI(address_gaussian,:).*FC(address_gaussian,:).*plx(address_gaussian,:),2)./sum(FC(address_gaussian,:).*plx(address_gaussian,:),2);
    % 对应rayleigh分布
    Mu(address_rayleigh, i+1) = sum(XI(address_rayleigh,:).^2.*FC(address_rayleigh,:).*plx(address_rayleigh,:),2)./(2*sum(FC(address_rayleigh,:).*plx(address_rayleigh,:),2));
    % 对应exponent分布
    Mu(address_exponent,i+1) = sum(FC(address_exponent,:).*plx(address_exponent,:),2)./(sum(XI(address_exponent,:).*FC(address_exponent,:).*plx(address_exponent,:),2));
    % 对应gaussian分布方差
    MU = M_Extand(Mu(:,i+1),K,L);
    Var(address_gaussian,i+1) = sum((XI(address_gaussian,:)-MU(address_gaussian,:)).^2.*FC(address_gaussian,:).*plx(address_gaussian,:),2)./sum(FC(address_gaussian,:).*plx(address_gaussian,:),2);
    
    % 固定分布不变
    W(address_fix, i+1) = W(address_fix, i);
    Mu(address_fix, i+1) = Mu(address_fix, i);
    Var(address_fix, i+1) = Var(address_fix, i);
    
    Error(i) = sum(abs(sum(pxl,1)-N_pixels/N));
end
VBN_EM = [Mu(:,upNum+1) sqrt(Var(:,upNum+1)) W(:,upNum+1)];
% [minError,i_num] = min(Error);
% VBN_EM = [Mu(:,i_num+1) sqrt(Var(:,i_num+1)) W(:,i_num+1)];

if Flag==1
    flag = {':b';'--c';'-.m';'-g';'-r';'-k';'-.k';'-.y';'--k';':k';':g'};
    
    EM_mu = zeros(K,1);
    EM_var = zeros(K,1);
    EM_sigma = zeros(K,1);
    Omega = zeros(K,1);
    RG_EM_curl = zeros(K,numel(N_pixels));
    plot(X_pixels, N_pixels/LengthRegion,'-k','LineWidth',1.5);xlim([0 200]); ylim([0, 0.1]);% 显示直方图曲线
    % axis([0 450 0 max(hc_Iout_new)]);
    grid on;hold on;
    for i = 1:K % 计算各类均值K_mean、均方差K_sigma百分比K_percent
        EM_mu(i) = VBN_EM(i,1);
        Omega(i) = VBN_EM(i,3);% 各子类分布曲线的最大值
        EM_var(i) =  VBN_EM(i,2)^2+eps(1);
        EM_sigma(i) = VBN_EM(i,2)+eps(1);
%         RG_EM_curl(i,:) = (i~=3)*Omega(i)*(EM_mu(i,1)).*exp(-(X_liver_vesselness.*EM_mu(i,1)))+...
%                        (i==3)*Omega(i)*(1/sqrt(2*pi)/EM_sigma(i)).*exp(-(X_liver_vesselness-EM_mu(i)).^2/(2*EM_var(i))); 
        RG_EM_curl(i,:) = (ismember(i, address_exponent))*Omega(i)*EM_mu(i).*exp(-(X_pixels.*EM_mu(i))) ...
                       + (ismember(i, address_gaussian))*Omega(i)*(1/sqrt(2*pi)/EM_sigma(i)).*exp(-(X_pixels-EM_mu(i)).^2/(2*EM_var(i))) ...
                       + (ismember(i, address_rayleigh))*Omega(i)*(X_pixels./EM_mu(i)^2).*exp(-(X_pixels.^2./(2*EM_mu(i)^2)));
        plot(X_pixels,RG_EM_curl(i,:),char(flag(i)),'LineWidth',1);
    end
    t = toc; disp(['using ' num2str(t) '秒']);
    legend_char = cell(K+2,1);
    legend_char{1} = char('Histogram of the filtering response data');
    for i = 1:K % % 编辑legend
%         if i<3
%            legend_char{1+i} = char(['Exponential ' num2str(i) ': lamit=' num2str((EM_mu(i)))...
%                ' w=' num2str(Omega(i))]);
%         else
%            legend_char{1+i} = char(['Gaussian ' num2str(i) ': mu=' num2str((EM_mu(i)))...
%                ' sigma=' num2str((EM_sigma(i))) ' w=' num2str(Omega(i))]);
%         end
        if ismember(i, address_exponent)
            legend_char{1+i} = char(['Exponential curl-line' num2str(i) ': lamit=' num2str((EM_mu(i)))...
              ' w=' num2str(Omega(i))]);
        elseif ismember(i, address_gaussian)
            legend_char{1+i} = char(['Gaussian curl-line' num2str(i) ': mu=' num2str(uint16(EM_mu(i)))...
              ' sigma=' num2str(uint16(EM_sigma(i))) ' w=' num2str(Omega(i))]);
        else
            legend_char{1+i} = char(['Rayleigh curl-line' num2str(i) ': sigma=' num2str((EM_mu(i)))...
              ' w=' num2str(Omega(i))]);
        end
    end
    plot(X_pixels,sum(RG_EM_curl,1),'-r','LineWidth',2);% 显示拟合后的曲线
    legend_char{K+2} = char('Fitting histogram');
    legend(legend_char{1:K+2});
    xlabel('Intensity');
    ylabel('Frequency');
    hold off
end
SumError = Error(dL);
end

%***********************************************************************************
%************** 子函数：求后验概率,f(k|xi) = wk*f(xi|k)/∑j=1:K(wj*f(xi|j))*********
function [plx,pxl] = pLX(xi,K,Mu,Var,W,address_gaussian,address_rayleigh,address_exponent)
% 计算后验概率矩阵plx
pxl = zeros(K,length(xi));% 初始化条件概率矩阵
plx = zeros(K,length(xi));% 初始化后验概率矩阵
% Var(1) = eps(1);  % 使得第一行方差为不等于零的最小数，以免1/sqrt(2*pi*Var(k))为NaN
% Var(2) = eps(1);
% pxl(1,:) =W(1)*(Mu(1)).*exp(-(xi.*Mu(1)));
% pxl(2,:) =W(2)*(Mu(2)).*exp(-(xi.*Mu(2)));
% pxl(3,:) =W(3)*(1/sqrt(2*pi*Var(3)))*exp(-(xi-Mu(3)).^2./(2*Var(3)));
Var(address_rayleigh) = eps(1);
Var(address_exponent) = eps(1);
pxl(address_rayleigh,:) = W(address_rayleigh).*(xi./Mu(address_rayleigh)).*exp(-(xi.^2./(2*Mu(address_rayleigh))));% 瑞利分布的权
%                                              ↓这个点超级重要！如果没有这个点会自动变成矩阵除法！天知道老子调了多久？！
pxl(address_gaussian,:) =W(address_gaussian).*(1./sqrt(2*pi*Var(address_gaussian))).*exp(-(xi-Mu(address_gaussian)).^2./(2*Var(address_gaussian))); % 高斯分布的权
pxl(address_exponent,:) =W(address_exponent).*(Mu(address_exponent)).*exp(-(xi.*Mu(address_exponent))); % 指数分布的权
% for k = 1:K               % 逐类别计算条件概率矩阵
%     pxl(k,:) = (k~=3)*W(k)*(Mu(k)).*exp(-(xi.*Mu(k))) + ...
%                (k==3)*W(k)*(1/sqrt(2*pi*Var(k)))*exp(-(xi-Mu(k)).^2./(2*Var(k)));
% end
Sum_pxl = sum(pxl,1)+eps(1);
for k = 1:K
    plx(k,:) = pxl(k,:)./Sum_pxl;
end
end

function [D] = M_Extand(Vector,K,L)
% size(Vector) = [K,1] or [1,L]
% 将Vector扩展为矩阵D，size(D)=[K,L]
D = zeros(K,L);
[a,b] = size(Vector);
if a>1 && b==1 %如果Vector是一列矢量K×1
   for j = 1:L
       D(:,j) = Vector;
   end    
end
if a==1 && b>1 %如果Vector是一行矢量1×L
    for i = 1:K
        D(i,:) = Vector;
    end
end
end