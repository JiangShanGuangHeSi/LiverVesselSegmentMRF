function [distribution, Object, VBN_Init, VBN_EM, Dx_MLE, sort_D, image_Dp] ...
    = DistributionInitialize_EM_MLE(distribution, kmeansStart, image, label_mask, label_fg_init, upNum, MleType, fig_visible)
% �Զ��ֲ���ʼ��-EM�㷨�Ż�-��������ʼ��� ȫ���̺���
% EM_kmeans_initialize
% EM_FMM
% ML_estimation
% ���룺
% distribution: ��ʾ�ֲ����ַ����飬���磬���������˹��distribution='ggg'���������ָ����һ��˹��distribution='eeg'
% kmeansStart: kmeans��ʼ���ֲ���label_fg_init = []ʱ����Ӧ����distribution��ͬ
% image: ͼ������������0-1֮���doubleͼ��Ҳ�������뷶Χ������int16ͼ��
% label_mask: ���д�����ͼ����ģ֮��
% label_fg_init: ���ֲַ�������ǰ�����֣����Ҫ����ǰ�������[]
% upNum: EM�㷨��������
% fig_visible: �Ƿ���ʾͼ��
if nargin==7
    fig_visible = false; % Ĭ�ϲ���ʾͼ��
end

Object = numel(distribution);

%% 1. ����Ԥ����ȡ������ֵ���������������ʵķ�Χ��������ֱ��ͼ
pixels_mask = image(label_mask);
% �����������ֵ��0-1�����ڵģ���������һ����������
if(all(pixels_mask > -1) && all(pixels_mask < 2)) % -1,2��Ϊ�˱���ʧ��
    Imax = 1000;
    pixels_mask = Imax * pixels_mask;
end
% ����ֲ�����e�ֲ��������ֲ���Ҫȥ���ֲ��е�0����������������
if (contains(distribution, 'e') || contains(distribution, 'r'))
    idx_pixels_mask_notzero = find(pixels_mask ~= 0);
    pixels_mask = pixels_mask(idx_pixels_mask_notzero);
end

%% 2. ��kmeans��ʼ��
figure(7)
fig7 = gcf;
if fig_visible == false
    set(fig7, 'Visible', 'off'); % ������
end
set(fig7,'Units','centimeter','Position',[5 5 25 10]);% ���ô��ڴ�С
subplot(1,2,1);
if isempty(label_fg_init) % ���û�б�������ʼ��
    VBN_Init = EM_kmeans_initialize(pixels_mask, distribution, kmeansStart);
else % ��������г�ʼ��
    hist_seg = GetDecomposedHist(image, label_fg_init, label_mask);
    VBN_Init_fg = [hist_seg.fg_mean, hist_seg.fg_std, double(sum(hist_seg.N_fg)) / sum(hist_seg.N_mask)];
    VBN_Init_bg = EM_kmeans_initialize(pixels_mask, distribution(1:end-1), kmeansStart); % ֻ��ʼ��������4����˹
    VBN_Init = [VBN_Init_bg; VBN_Init_fg];
end

%% 3. EM�㷨���Ʒֲ�
figure(7);
if fig_visible == false
    set(fig7, 'Visible', 'off'); % ������
end
subplot(1,2,2);
if isempty(label_fg_init) % ���û�б�������ʼ��
    [VBN_EM, SumError] = EM_FMM(pixels_mask,VBN_Init,upNum,distribution,1);
else % ��������г�ʼ��
    [VBN_EM, SumError] = EM_FMM(pixels_mask,VBN_Init,upNum,distribution,1,Object);% �������ֲ�����
end

%% 4. ���������
pixels_mask = image(label_mask); % ֮ǰ�Ѿ��޸Ĺ������Ҫ�ָ���Χ
if(all(pixels_mask > -1) && all(pixels_mask < 2)) % -1,2��Ϊ�˱���ʧ��
    [Dx_MLE, sort_D, image_Dp] = ML_estimation(image * Imax,VBN_EM,Object,label_mask,distribution, MleType);  % ������Ȼ���Ƶĳ�ʼ��ǳ�
else
    [Dx_MLE, sort_D, image_Dp] = ML_estimation(image,VBN_EM,Object,label_mask,distribution, MleType);
end

end