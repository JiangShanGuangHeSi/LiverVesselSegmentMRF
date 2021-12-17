function [label_large,label_small] = SeparationVesselSize(label,r)
%VesselSizeSeparation ����Ѫ��label�еĴ�СѪ��
%   �ο���An Attention-guided Deep Neural Network with Multi-scale Feature Fusion for Liver Vessel Segmentation
%   �����ǣ��ô�СΪd�ĺ˽��п����㣬ֱ��<=d-1��Ѫ�ܻᱻȥ��
SE = strel('sphere', r);
label_large = imopen(label, SE);

% % ���Ľ���
% SE2 = strel('sphere', r+2); % ����������
% label_large_dilate = imdilate(label_large, SE2);
% ��ԭ�桿
label_large_dilate = label_large;

label_large = and(logical(label), logical(label_large_dilate));
label_small = and(logical(label), ~logical(label_large));
end