function [box] = GetBox3d(image)
%GetBox ������άͼ��Ϊ0������Χ
%   ����ֵ��
%       box ��i��Ϊ��i��ά�ȵķ�Χ��ʼֵ����ֵֹ
box = zeros(3, 2);
vertical_indicies = find(any(image, [2, 3])); % ��һ��ά�ȷ�������
horizontal_indicies = find(any(image, [1, 3])); % �ڶ���ά�ȷ�������
slices_indicies = find(any(image, [1, 2])); % ������ά�ȷ�������
box(1, 1) = vertical_indicies(1);
box(1, 2) = vertical_indicies(end);
box(2, 1) = horizontal_indicies(1);
box(2, 2) = horizontal_indicies(end);
box(3, 1) = slices_indicies(1);
box(3, 2) = slices_indicies(end);
end

