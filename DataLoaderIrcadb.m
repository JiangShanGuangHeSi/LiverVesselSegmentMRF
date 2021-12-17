function [image_mask, image_mask_norm, label_mask, image_box, spacing, slice_show, label_save, label_fg] ...
    = DataLoaderIrcadb(dict_path)
% ���ܣ���ȡircadb nii��ʽ�����ݼ����ü�������
% image_mask, image_mask_norm
% label_fg, label_bg, label_mask
% pixel_mask, pixel_fg, pixel_bg
% idx_mask, idx_fg, idx_bg
% ����image_box������spacing����ʾ�õ�slice_show���洢�õ�label_save

% ����������޸�Ϊdict�ṹ��������ֵlabel_fg�ɵڶ����ƶ������һ����2021.8��

dataset_root = dict_path.root_dataset;
image_path = dict_path.path_image;
label_path_liver = dict_path.path_label_liver;
label_path_port = dict_path.path_label_port;
label_path_vena = dict_path.path_label_vena;


%% ����ͼ��
image = load_nii(image_path);
spacing = image.hdr.dime.pixdim(2:4); % ����spacing
image = image.img;
label_liver = load_nii(label_path_liver);
label_save = label_liver; % ��һ��������
label_liver = label_liver.img;
label_port = load_nii(label_path_port);
label_port = label_port.img;
label_vena = load_nii(label_path_vena);
label_vena = label_vena.img;

%% ����box
image_box = GetBox3d(label_liver); % boxҪ��ȥ����֮ǰ�㣬��Ȼ������

%% ������ȥ����������Щ�����У���case17�������ǳ�Ӱ����
% ������ȥ��case21������ΪHyperplasie������
label_tumor = false(size(image));
case_file_list = dir(dataset_root);
for i = 1:length(case_file_list)
    if~case_file_list(i).isdir && (contains(case_file_list(i).name, 'tumor') || contains(case_file_list(i).name, 'Hyperplasie'))
        label_tumor_temp = load_nii([dataset_root, case_file_list(i).name]);
        label_tumor = label_tumor | logical(label_tumor_temp.img);
    end
end
%     figure(1); imshow(max(label_tumor, [], 3));
label_liver = uint8(logical(label_liver) & ~label_tumor);

image = single(image);

%% �ϲ��ž����͸ξ�����ǩ�������ڸ���mask��
label_fg = label_liver ~= 0 & (label_port ~= 0 | label_vena ~= 0);% �����ڲ�Ѫ�ܱ�ǩ����
label_bg = label_liver ~= 0 & ~(label_port ~= 0 | label_vena ~= 0);% �����ڲ���������Ѫ�ܣ���ǩ����

%% �ü�
label_fg = label_fg(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2));
label_bg = label_bg(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2));
label_mask = logical(label_liver(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2)));
image_mask = image(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2));
image_mask(~label_mask) = 0;

%% ͼ��Ӵ���40 350�Ǹ�ǻ�����Ĵ���ֻ��Ϊ����ʾ�ã�
WINDOW_POSITION = 40;
WINDOW_WIDTH = 350;
Imax = WINDOW_POSITION + WINDOW_WIDTH / 2;
Imin = WINDOW_POSITION - WINDOW_WIDTH / 2;
image_mask_norm = (image_mask - Imin) / (Imax - Imin);

%% ������ʾ�Ĳ���
slice_show = round((image_box(3,2) - image_box(3,1)) / 2);

end