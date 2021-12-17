function [label_large,label_small] = SeparationVesselSize(label,r)
%VesselSizeSeparation 分离血管label中的大小血管
%   参考：An Attention-guided Deep Neural Network with Multi-scale Feature Fusion for Liver Vessel Segmentation
%   大意是：用大小为d的核进行开运算，直径<=d-1的血管会被去除
SE = strel('sphere', r);
label_large = imopen(label, SE);

% % 【改进】
% SE2 = strel('sphere', r+2); % 再膨胀两下
% label_large_dilate = imdilate(label_large, SE2);
% 【原版】
label_large_dilate = label_large;

label_large = and(logical(label), logical(label_large_dilate));
label_small = and(logical(label), ~logical(label_large));
end