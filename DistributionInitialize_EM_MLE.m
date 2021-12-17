function [distribution, Object, VBN_Init, VBN_EM, Dx_MLE, sort_D, image_Dp] ...
    = DistributionInitialize_EM_MLE(distribution, kmeansStart, image, label_mask, label_fg_init, upNum, MleType, fig_visible)
% 自动分布初始化-EM算法优化-最大后验概率计算 全流程函数
% EM_kmeans_initialize
% EM_FMM
% ML_estimation
% 输入：
% distribution: 表示分布的字符数组，比如，如果是三高斯，distribution='ggg'，如果是两指数和一高斯，distribution='eeg'
% kmeansStart: kmeans初始化分布，label_fg_init = []时长度应该与distribution相同
% image: 图像本身，可以输入0-1之间的double图像，也可以输入范围正常的int16图像
% label_mask: 所有处理在图像掩模之内
% label_fg_init: 保持分布不动的前景部分，如果要估计前景输入空[]
% upNum: EM算法迭代次数
% fig_visible: 是否显示图像
if nargin==7
    fig_visible = false; % 默认不显示图像
end

Object = numel(distribution);

%% 1. 数据预处理，取出像素值，将它放缩到合适的范围，并计算直方图
pixels_mask = image(label_mask);
% 如果所有像素值是0-1区间内的，对它进行一个坐标拉伸
if(all(pixels_mask > -1) && all(pixels_mask < 2)) % -1,2是为了避免失误
    Imax = 1000;
    pixels_mask = Imax * pixels_mask;
end
% 如果分布包含e分布或瑞利分布，要去除分布中的0分量，避免计算错误
if (contains(distribution, 'e') || contains(distribution, 'r'))
    idx_pixels_mask_notzero = find(pixels_mask ~= 0);
    pixels_mask = pixels_mask(idx_pixels_mask_notzero);
end

%% 2. 用kmeans初始化
figure(7)
fig7 = gcf;
if fig_visible == false
    set(fig7, 'Visible', 'off'); % 不弹出
end
set(fig7,'Units','centimeter','Position',[5 5 25 10]);% 设置窗口大小
subplot(1,2,1);
if isempty(label_fg_init) % 如果没有背景景初始化
    VBN_Init = EM_kmeans_initialize(pixels_mask, distribution, kmeansStart);
else % 如果背景有初始化
    hist_seg = GetDecomposedHist(image, label_fg_init, label_mask);
    VBN_Init_fg = [hist_seg.fg_mean, hist_seg.fg_std, double(sum(hist_seg.N_fg)) / sum(hist_seg.N_mask)];
    VBN_Init_bg = EM_kmeans_initialize(pixels_mask, distribution(1:end-1), kmeansStart); % 只初始化背景的4个高斯
    VBN_Init = [VBN_Init_bg; VBN_Init_fg];
end

%% 3. EM算法估计分布
figure(7);
if fig_visible == false
    set(fig7, 'Visible', 'off'); % 不弹出
end
subplot(1,2,2);
if isempty(label_fg_init) % 如果没有背景景初始化
    [VBN_EM, SumError] = EM_FMM(pixels_mask,VBN_Init,upNum,distribution,1);
else % 如果背景有初始化
    [VBN_EM, SumError] = EM_FMM(pixels_mask,VBN_Init,upNum,distribution,1,Object);% 保持最后分布不动
end

%% 4. 最大后验估计
pixels_mask = image(label_mask); % 之前已经修改过，因此要恢复范围
if(all(pixels_mask > -1) && all(pixels_mask < 2)) % -1,2是为了避免失误
    [Dx_MLE, sort_D, image_Dp] = ML_estimation(image * Imax,VBN_EM,Object,label_mask,distribution, MleType);  % 极大似然估计的初始标记场
else
    [Dx_MLE, sort_D, image_Dp] = ML_estimation(image,VBN_EM,Object,label_mask,distribution, MleType);
end

end