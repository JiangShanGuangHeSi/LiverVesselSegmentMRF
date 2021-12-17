function [VBN_Init] = EM_kmeans_initialize(pixels, distribution, kmeansStart)
% ���ݷֲ�����distribution����ͼ������ֵpixels�ֲ�������ʼ��
% ע����������pixels��Ҫ����ɢֵ��������С������������ֲ��а���e�ֲ�����Ҫע������ֵӦ���ϸ����0��������
% pixels�е����������kmeans����

% �����������
K = numel(distribution);
% ����distribution�ַ�������ֲ���λ��˳��
address_gaussian = find(distribution == 'g');
address_rayleigh = find(distribution == 'r');
address_exponent = find(distribution == 'e');

% ����ֱ��ͼ
% �����e�ֲ���һ��Ҫ��1��ʼ����1��ʼ����1��ʼ��
[N, X] = hist(pixels, min(pixels):max(pixels)); 
 
% [idx,ctrs] = kmeans(pixels(:),K,'start',[a,b,c]);
if nargin==2 || (nargin==3 && isempty(kmeansStart))
    [idx,ctrs] = kmeans(pixels(:),K);
elseif nargin==3 && ~isempty(kmeansStart)
    [idx,ctrs] = kmeans(pixels(:),K,'Start',kmeansStart);
end
[idx,ctrs] = Kmean_reorder(idx,ctrs);% ���ջҶ��������ɵ����ߵ�˳���������idx��ctrs

tic;
K_mu = zeros(K, 1);
K_var = zeros(K,1);
K_sigma = zeros(K,1);
Omega = zeros(K,1);
RG_K_curl = zeros(K,numel(N));
% LengthRegion = length(liver_idx);
% LengthRegion = length(pixels_liver_vesselness_notzero);
LengthRegion = numel(pixels);
plot(X,N/LengthRegion,'-k', 'LineWidth',2);xlim([0 200]); ylim([0, 0.1]);% ��ʾֱ��ͼ����
% subplot(1,2,2);plot(X_liver_vesselness,N_liver_vesselness/LengthRegion,'-k', 'LineWidth',2);% ��ʾֱ��ͼ����
% axis([0 450 0 max(hc_Iout_new)]);
grid on;hold on;
flag = {':b';'--c';'-.m';'-g';'-r';'-k';'-.k';'-.y';'--k';':k';':g'};
for i = 1:K % ������������K_var���ٷֱ�K_percent�����˹����ͼ
    % 1. ����Ȩ�أ�����ں���д��w
    Omega(i) = length(find(idx==i))/LengthRegion;% ������ֲ����ߵ����ֵ
    % 2. ����mu
%         K_mu(i) = (i~=K)*1.0/mean(pixels_liver_vesselness_notzero(idx==i)) + (i==K)*mean(pixels_liver_vesselness_notzero(idx==i));
    % ˵��һ�£���˹�Ħ̵��ڷֲ����������������ֲ����������ڦҡ�(��/2)��ָ���ֲ��Ħ�Ϊ������֮һ
    K_mu(i) = (ismember(i, address_exponent))*1.0/mean(pixels(idx==i)) ...
            + (ismember(i, address_gaussian))*mean(pixels(idx==i)) ...
            + (ismember(i, address_rayleigh)*mean(pixels(idx==i)))*sqrt(2/pi); % ��������ڶ�����������ǣ������ֲ����ܲ�����ô��ʼ���ģ�
    % 3. ����sigma��Ŀǰֻ��Gaussian�ֲ���������������壩
    K_var(i) = var(pixels(idx==i));
    K_sigma(i) = sqrt(K_var(i));
    % 4. ��������
%         RG_K_curl(i,:) = (i~=K)*Omega(i)*K_mu(i).*exp(-(X_liver_vesselness.*K_mu(i)))+...
%                        (i==K)*Omega(i)*(1/sqrt(2*pi)/K_sigma(i)).*exp(-(X_liver_vesselness-K_mu(i)).^2/(2*K_var(i)));
    RG_K_curl(i,:) = (ismember(i, address_exponent))*Omega(i)*K_mu(i).*exp(-(X.*K_mu(i))) ...
                   + (ismember(i, address_gaussian))*Omega(i)*(1/sqrt(2*pi)/K_sigma(i)).*exp(-(X-K_mu(i)).^2/(2*K_var(i))) ...
                   + (ismember(i, address_rayleigh))*Omega(i)*(X./K_mu(i)^2).*exp(-(X.^2./(2*K_mu(i)^2)));
    plot(X,RG_K_curl(i,:),char(flag(i)),'LineWidth',1);%���Ƹ�����ֲ�����
end
t = toc; disp(['using ' num2str(t) '��']);


legend_char = cell(K+2,1);
legend_char{1} = char('Original histogram');
for i = 1:K % �༭legend
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
plot(X,sum(RG_K_curl,1),'-r','LineWidth',2);% ��ʾ��Ϻ������
legend_char{K+2} = char('Init-fitting histogram');
legend(legend_char{1:K+2});
xlabel('Intensity');
ylabel('Frequency');
hold off

VBN_Init = [K_mu, K_sigma, Omega];