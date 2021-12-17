function [decomposed_hist] = GetDecomposedHist(image, label, mask)
% 绘制图像mask范围内的直方图，label标记的前景的直方图和背景的直方图
% 返回图片的句柄
decomposed_hist = {};

idx_mask = mask ~= 0; % 分析范围的像素索引
pixels_mask = image(idx_mask); % 分析范围的像素
if ~isempty(label)
    idx_fg = mask ~=0 & label ~=0; % 前景索引
    pixels_fg = image(idx_fg); % 前景像素
    idx_bg = mask ~= 0 & label == 0; % 背景索引
    pixels_bg = image(idx_bg); % 背景像素
end

% 计算直方图
[N_mask, X_mask] = hist(pixels_mask, min(pixels_mask):max(pixels_mask));
decomposed_hist.N_mask = N_mask;
decomposed_hist.X_mask = X_mask;
if ~isempty(label)
    [N_fg, X_fg] = hist(pixels_fg, min(pixels_mask):max(pixels_mask));
    decomposed_hist.N_fg = N_fg;
    decomposed_hist.X_fg = X_fg;
    [N_bg, X_bg] = hist(pixels_bg, min(pixels_mask):max(pixels_mask));
    decomposed_hist.N_bg = N_bg;
    decomposed_hist.X_bg = X_bg;
end

%分析直方图
decomposed_hist.mean = mean(pixels_mask);
decomposed_hist.std = std(pixels_mask);
if ~isempty(label)
    decomposed_hist.fg_mean = mean(pixels_fg);
    decomposed_hist.fg_std = std(pixels_fg);
    decomposed_hist.bg_mean = mean(pixels_bg);
    decomposed_hist.bg_std = std(pixels_bg);
    % 统计前景均值在整个直方图的分位数（前百分之几）
    quantile_fg_mean = sum(pixels_mask > mean(pixels_fg)) / numel(pixels_mask);
    decomposed_hist.quantile_fg_mean = quantile_fg_mean;
end
% % 显示
% figure(fig_id); hold on
% % 设置窗口大小
% h_fig = gcf;
% set(gcf,'Units','centimeter','Position',[5 5 10 7]);
% plot(X_mask, N_mask); 
% plot(X_fg, N_fg);
% plot(X_bg, N_bg);
% xlim([-200 400]); ylim([0, 7e4]);
% legend('liver', 'vessel', 'background'); hold off;
%     saveas(fig1, [dir_name, 'figure5, hist of origin image'], 'fig');
% saveas(h_fig, [dir_name, 'figure6, histogram specification'], 'fig');
end