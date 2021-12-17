function [Dout,sort_D, Dp] = ML_estimation(A,VBN,Object,flagIMG,distribution,type)
% �������壺���м��������Ȼ����
% A,����������Ѫ�ܺͱ���������ݣ�
% VBN3�������ֵ�ͷ��
% W��������ϵ��
% sort_D:��Ȼ���ʺ�ͷ������IndexA�ĺϳɾ�����length(IndexA)�У�Object+1��
tic;
disp('ML_estimation ...');

% ����distribution�ַ�������ֲ���λ��˳��
address_gaussian = find(distribution == 'g');
address_rayleigh = find(distribution == 'r');
address_exponent = find(distribution == 'e');

%��20210907����,�Ľ��˸���������⡿
AM = A(A>=(VBN(Object,1) + VBN(Object,2)));
A(A>=(VBN(Object,1) + VBN(Object,2))) = min(AM(:));
Am = A(A<=(VBN(1,1) - VBN(1,2)));
A(A<=(VBN(1,1) - VBN(1,2))) = max(Am(:));
clear AM Am;

Dout = zeros(size(A));
Dp = Dout;
IndexA = find(flagIMG>=1);%ȡ��flagIMG���Ϊ1����������
% IndexA = Index_head;
N = numel(IndexA);
L = size(VBN,1);
A = repmat(A(IndexA)',L,1);%A���ΪL��1�е�A(:)'����
mu = repmat(VBN(:,1),1,N);%����L��Ni�е�mu����
sigma = repmat(VBN(:,2),1,N);%����L��Ni�е�sigma����
W = repmat(VBN(:,3),1,N);%����L��Ni�е�W����
Li = find(1:L~=Object)';% ����1��K-1��

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
        DMPL_vector = ((pxl_k(Object,:)./W(Object,:))>(sum(pxl_k(Li,:),1)./sum(W(Li,:),1)));%��һ�֣�Ȩƽ��Լ��Ŀ��ͱ�����
        Dp(IndexA) = (pxl_k(Object,:)./W(Object,:)) ./ (pxl_k(Object,:)./W(Object,:) + sum(pxl_k(Li,:),1)./sum(W(Li,:),1)); % ���������ĸ���ͼ
    case 2
        DMPL_vector = Object*(pxl_k(Object,:)>max(pxl_k(Li,:),[],1));%�ڶ��֣�Ѫ������ڱ���������ֵ
        Dp(IndexA) = pxl_k(Object,:)./(pxl_k(Object,:)+max(pxl_k(Li,:),[],1));
    case 3
        DMPL_vector = Object*(pxl_k(Object,:)>sum(pxl_k(Li,:),1)/(L-1));%�����֣�Ѫ������ڱ�����ľ�ֵ,��������һ���㷨Ч������
        Dp(IndexA) = pxl_k(Object,:)./(sum(pxl_k,1));
end
Dout(IndexA)= DMPL_vector';%��ʸ�����������ֵ
index_PLX_object = [pxl_k' IndexA];% ���flagIMG==1�Ŀռ����أ������������飬��һ����pxl_k(Object)���ڶ��������صı��
sort_D = sortrows(index_PLX_object,-3);%��ȥ�ԿǺ�Ŀռ��У�����Ŀ���ࣨpxl_k��:,4�����ɸ����͵�˳���������
t = toc;
disp(['finished, time consumption of this step is ' num2str(t) ' s']);