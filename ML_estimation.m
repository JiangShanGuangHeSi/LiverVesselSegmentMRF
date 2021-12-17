function [Dout,sort_D, Dp] = ML_estimation(A,VBN,Object,flagIMG,distribution,type)
% 函数意义：并行计算最大似然估计
% A,噪声条件下血管和背景混合数据；
% VBN3：各类均值和方差；
% W：类别比例系数
% sort_D:似然概率和头部索引IndexA的合成矩阵，有length(IndexA)行，Object+1列
tic;
disp('ML_estimation ...');

% 根据distribution字符串计算分布的位置顺序
address_gaussian = find(distribution == 'g');
address_rayleigh = find(distribution == 'r');
address_exponent = find(distribution == 'e');

%【20210907修正,改进了概率溢出问题】
AM = A(A>=(VBN(Object,1) + VBN(Object,2)));
A(A>=(VBN(Object,1) + VBN(Object,2))) = min(AM(:));
Am = A(A<=(VBN(1,1) - VBN(1,2)));
A(A<=(VBN(1,1) - VBN(1,2))) = max(Am(:));
clear AM Am;

Dout = zeros(size(A));
Dp = Dout;
IndexA = find(flagIMG>=1);%取被flagIMG标记为1的像素索引
% IndexA = Index_head;
N = numel(IndexA);
L = size(VBN,1);
A = repmat(A(IndexA)',L,1);%A阵变为L行1列的A(:)'矩阵
mu = repmat(VBN(:,1),1,N);%产生L行Ni列的mu矩阵
sigma = repmat(VBN(:,2),1,N);%产生L行Ni列的sigma矩阵
W = repmat(VBN(:,3),1,N);%产生L行Ni列的W矩阵
Li = find(1:L~=Object)';% 产生1至K-1行

% pxl_k = [W(1:2,:).*(mu(1:2,:)).*(exp(-(A(1:2,:).*mu(1:2,:))));...
%          W(3,:).*(1./sqrt(2*pi*sigma(3,:).^2)).*exp(-(A(3,:)-mu(3,:)).^2./(2*sigma(3,:).^2))];
pxl_k = zeros(L, N);
pxl_k(address_gaussian, :) = W(address_gaussian,:).*(1./sqrt(2*pi*sigma(address_gaussian,:).^2)).*exp(-(A(address_gaussian,:)-mu(address_gaussian,:)).^2./(2*sigma(address_gaussian,:).^2));
pxl_k(address_exponent, :) = W(address_exponent,:).*(mu(address_exponent,:)).*(exp(-(A(address_exponent,:).*mu(address_exponent,:))));
pxl_k(address_rayleigh, :) = W(address_rayleigh,:).*((A(address_rayleigh,:)./mu(address_rayleigh,:))).*(exp(-(A(address_rayleigh,:).^2./(2*mu(address_rayleigh,:)))));

if nargin == 5
    type = 1; % 
end
switch type
    case 1
        DMPL_vector = ((pxl_k(Object,:)./W(Object,:))>(sum(pxl_k(Li,:),1)./sum(W(Li,:),1)));%第一种：权平均约束目标和背景类
        Dp(IndexA) = (pxl_k(Object,:)./W(Object,:)) ./ (pxl_k(Object,:)./W(Object,:) + sum(pxl_k(Li,:),1)./sum(W(Li,:),1)); % 返回真正的概率图
    case 2
        DMPL_vector = Object*(pxl_k(Object,:)>max(pxl_k(Li,:),[],1));%第二种：血管类大于背景类的最大值
        Dp(IndexA) = pxl_k(Object,:)./(pxl_k(Object,:)+max(pxl_k(Li,:),[],1));
    case 3
        DMPL_vector = Object*(pxl_k(Object,:)>sum(pxl_k(Li,:),1)/(L-1));%第三种：血管类大于背景类的均值,和上述第一种算法效果类似
        Dp(IndexA) = pxl_k(Object,:)./(sum(pxl_k,1));
end
Dout(IndexA)= DMPL_vector';%列矢量向输出矩阵赋值
index_PLX_object = [pxl_k' IndexA];% 针对flagIMG==1的空间体素，构造两列数组，第一列是pxl_k(Object)，第二列是体素的编号
sort_D = sortrows(index_PLX_object,-3);%在去脑壳后的空间中，按照目标类（pxl_k（:,4））由高至低的顺序输出矩阵
t = toc;
disp(['finished, time consumption of this step is ' num2str(t) ' s']);