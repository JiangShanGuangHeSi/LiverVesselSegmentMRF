function [VBN_Init] = EM_kmeans_initialize(pixels, distribution, kmeansStart)
% 根据分布类型distribution进行图像像素值pixels分布参数初始化
% 注意输入像素pixels需要是离散值，不能是小数！！！如果分布中包含e分布，需要注意像素值应该严格大于0！！！。
% pixels中的像素类别用kmeans聚类

% 计算类别数量
K = numel(distribution);
% 根据distribution字符串计算分布的位置顺序
address_gaussian = find(distribution == 'g');
address_rayleigh = find(distribution == 'r');
address_exponent = find(distribution == 'e');

% 计算直方图
% 如果有e分布，一定要从1开始，从1开始，从1开始！
[N, X] = hist(pixels, min(pixels):max(pixels)); 
 
% [idx,ctrs] = kmeans(pixels(:),K,'start',[a,b,c]);
if nargin==2 || (nargin==3 && isempty(kmeansStart))
    [idx,ctrs] = kmeans(pixels(:),K);
elseif nargin==3 && ~isempty(kmeansStart)
    [idx,ctrs] = kmeans(pixels(:),K,'Start',kmeansStart);
end
[idx,ctrs] = Kmean_reorder(idx,ctrs);% 按照灰度类中心由低至高的顺序重新输出idx和ctrs

tic;
K_mu = zeros(K, 1);
K_var = zeros(K,1);
K_sigma = zeros(K,1);
Omega = zeros(K,1);
RG_K_curl = zeros(K,numel(N));
% LengthRegion = length(liver_idx);
% LengthRegion = length(pixels_liver_vesselness_notzero);
LengthRegion = numel(pixels);
plot(X,N/LengthRegion,'-k', 'LineWidth',2);xlim([0 200]); ylim([0, 0.1]);% 显示直方图曲线
% subplot(1,2,2);plot(X_liver_vesselness,N_liver_vesselness/LengthRegion,'-k', 'LineWidth',2);% 显示直方图曲线
% axis([0 450 0 max(hc_Iout_new)]);
grid on;hold on;
flag = {':b';'--c';'-.m';'-g';'-r';'-k';'-.k';'-.y';'--k';':k';':g'};
for i = 1:K % 计算各类均方差K_var、百分比K_percent并绘高斯曲线图
    % 1. 计算权重，这个在后面写作w
    Omega(i) = length(find(idx==i))/LengthRegion;% 各子类分布曲线的最大值
    % 2. 计算mu
%         K_mu(i) = (i~=K)*1.0/mean(pixels_liver_vesselness_notzero(idx==i)) + (i==K)*mean(pixels_liver_vesselness_notzero(idx==i));
    % 说明一下，高斯的μ等于分布的期望，和瑞利分布的期望等于σ√(Π/2)，指数分布的λ为期望分之一
    K_mu(i) = (ismember(i, address_exponent))*1.0/mean(pixels(idx==i)) ...
            + (ismember(i, address_gaussian))*mean(pixels(idx==i)) ...
            + (ismember(i, address_rayleigh)*mean(pixels(idx==i)))*sqrt(2/pi); % 这句适用于多种情况，但是！瑞利分布可能不是这么初始化的！
    % 3. 计算sigma（目前只有Gaussian分布的这个参数有意义）
    K_var(i) = var(pixels(idx==i));
    K_sigma(i) = sqrt(K_var(i));
    % 4. 计算曲线
%         RG_K_curl(i,:) = (i~=K)*Omega(i)*K_mu(i).*exp(-(X_liver_vesselness.*K_mu(i)))+...
%                        (i==K)*Omega(i)*(1/sqrt(2*pi)/K_sigma(i)).*exp(-(X_liver_vesselness-K_mu(i)).^2/(2*K_var(i)));
    RG_K_curl(i,:) = (ismember(i, address_exponent))*Omega(i)*K_mu(i).*exp(-(X.*K_mu(i))) ...
                   + (ismember(i, address_gaussian))*Omega(i)*(1/sqrt(2*pi)/K_sigma(i)).*exp(-(X-K_mu(i)).^2/(2*K_var(i))) ...
                   + (ismember(i, address_rayleigh))*Omega(i)*(X./K_mu(i)^2).*exp(-(X.^2./(2*K_mu(i)^2)));
    plot(X,RG_K_curl(i,:),char(flag(i)),'LineWidth',1);%绘制各子类分布曲线
end
t = toc; disp(['using ' num2str(t) '秒']);


legend_char = cell(K+2,1);
legend_char{1} = char('Original histogram');
for i = 1:K % 编辑legend
    if ismember(i, address_exponent)
        legend_char{1+i} = char(['Exponential curl-line' num2str(i) ': lamit=' num2str((K_mu(i)))...
          ' w=' num2str(Omega(i))]);
    elseif ismember(i, address_gaussian)
        legend_char{1+i} = char(['Gaussian curl-line' num2str(i) ': mu=' num2str(uint16(K_mu(i)))...
          ' sigma=' num2str(uint16(K_sigma(i))) ' w=' num2str(Omega(i))]);
    else
        legend_char{1+i} = char(['Rayleigh curl-line' num2str(i) ': sigma=' num2str((K_mu(i)))...
          ' w=' num2str(Omega(i))]);
    end
end
plot(X,sum(RG_K_curl,1),'-r','LineWidth',2);% 显示拟合后的曲线
legend_char{K+2} = char('Init-fitting histogram');
legend(legend_char{1:K+2});
xlabel('Intensity');
ylabel('Frequency');
hold off

VBN_Init = [K_mu, K_sigma, Omega];