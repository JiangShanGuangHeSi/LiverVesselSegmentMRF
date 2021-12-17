function [image_mask, image_mask_norm, label_mask, image_box, spacing, slice_show, label_save, label_fg] ...
    = DataLoaderZDYY(dict_path)
% ���ܣ���ȡnii��ʽ�����ݼ����ü�������
% image_mask, image_mask_norm
% label_fg, label_bg, label_mask
% pixel_mask, pixel_fg, pixel_bg
% idx_mask, idx_fg, idx_bg
% ����image_box������spacing����ʾ�õ�slice_show���洢�õ�label_save

% ��dataloader2�Ļ����Ͻ������ĳ�һ��dict�ṹ��2021.8��

image_path = dict_path.path_data_nii;
label_path_liver = dict_path.path_mask_nii;
% if(exist('dict_path.path_label_nii', 'var') && ~isempty(dict_path.path_label_nii)) % �жϽṹ��ĳ���ֶ��Ƿ����Ӧ����isfield(dict_path, 'path_label_nii')
if(isfield(dict_path, 'path_label_nii') && ~isempty(dict_path.path_label_nii) && exist(dict_path.path_label_nii, 'file'))
    do_fg = true;
    label_path_vessel = dict_path.path_label_nii;
else
    do_fg = false;
end
idx_mask = dict_path.idx_mask;

%% ����ͼ��
image = load_nii(image_path);
spacing = image.hdr.dime.pixdim(2:4); % ����spacing
image = image.img;
label_liver = load_nii(label_path_liver);
label_save = label_liver; % ��һ��������
label_liver = label_liver.img;
if do_fg
    label_fg = load_nii(label_path_vessel);
    label_fg = label_fg.img;
else
    label_fg = [];
end

%% ��ATTENTION: ��Բ�ͬ���ݼ�Ҫ�ò�ͬ�ķ�ʽ��ȡlabel��
% �ؽ��������ݼ���liver��ע = 6
% nnunetԤ�ָ���liver��ע = 1
label_liver = label_liver==idx_mask;
% �ҵ����ݼ�����Щvessel��ע= 1����Щvessel��ע�ξ���Ϊ1���ž���Ϊ0
if do_fg
    label_fg = label_fg ~= 0;
    label_fg = logical(label_fg .* label_liver);
end

%% ����box
image_box = GetBox3d(label_liver); % boxҪ��ȥ����֮ǰ�㣬��Ȼ������

%% �ü�
image_mask = image(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2));
% image_mask(label_mask) = 0;
label_mask = logical(label_liver(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2)));
if do_fg
    label_fg = logical(label_fg(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2)));
end

%% ͼ��Ӵ���40 350�Ǹ�ǻ�����Ĵ���ֻ��Ϊ����ʾ�ã�
WINDOW_POSITION = 40;
WINDOW_WIDTH = 350;
Imax = WINDOW_POSITION + WINDOW_WIDTH / 2;
Imin = WINDOW_POSITION - WINDOW_WIDTH / 2;
image_mask_norm = (double(image_mask) - Imin) / (Imax - Imin);

%% ������ʾ�Ĳ���
slice_show = round((image_box(3,2) - image_box(3,1)) / 2);

end