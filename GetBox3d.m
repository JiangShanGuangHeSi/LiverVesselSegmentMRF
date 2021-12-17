function [box] = GetBox3d(image)
%GetBox 返回三维图像不为0的区域范围
%   返回值：
%       box 第i行为第i个维度的范围初始值和终止值
box = zeros(3, 2);
vertical_indicies = find(any(image, [2, 3])); % 第一个维度非零数组
horizontal_indicies = find(any(image, [1, 3])); % 第二个维度非零数组
slices_indicies = find(any(image, [1, 2])); % 第三个维度非零数组
box(1, 1) = vertical_indicies(1);
box(1, 2) = vertical_indicies(end);
box(2, 1) = horizontal_indicies(1);
box(2, 2) = horizontal_indicies(end);
box(3, 1) = slices_indicies(1);
box(3, 2) = slices_indicies(end);
end

