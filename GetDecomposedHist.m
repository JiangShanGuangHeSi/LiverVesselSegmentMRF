function [decomposed_hist] = GetDecomposedHist(image, label, mask)
% ����ͼ��mask��Χ�ڵ�ֱ��ͼ��label��ǵ�ǰ����ֱ��ͼ�ͱ�����ֱ��ͼ
% ����ͼƬ�ľ��
decomposed_hist = {};

idx_mask = mask ~= 0; % ������Χ����������
pixels_mask = image(idx_mask); % ������Χ������
if ~isempty(label)
    idx_fg = mask ~=0 & label ~=0; % ǰ������
    pixels_fg = image(idx_fg); % ǰ������
    idx_bg = mask ~= 0 & label == 0; % ��������
    pixels_bg = image(idx_bg); % ��������
end

% ����ֱ��ͼ
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

%����ֱ��ͼ
decomposed_hist.mean = mean(pixels_mask);
decomposed_hist.std = std(pixels_mask);
if ~isempty(label)
    decomposed_hist.fg_mean = mean(pixels_fg);
    decomposed_hist.fg_std = std(pixels_fg);
    decomposed_hist.bg_mean = mean(pixels_bg);
    decomposed_hist.bg_std = std(pixels_bg);
    % ͳ��ǰ����ֵ������ֱ��ͼ�ķ�λ����ǰ�ٷ�֮����
    quantile_fg_mean = sum(pixels_mask > mean(pixels_fg)) / numel(pixels_mask);
    decomposed_hist.quantile_fg_mean = quantile_fg_mean;
end
% % ��ʾ
% figure(fig_id); hold on
% % ���ô��ڴ�С
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