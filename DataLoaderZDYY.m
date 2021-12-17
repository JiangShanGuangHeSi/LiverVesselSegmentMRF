function [image_mask, image_mask_norm, label_mask, image_box, spacing, slice_show, label_save, label_fg] ...
    = DataLoaderZDYY(dict_path)
% 功能：读取nii格式的数据集，裁剪，返回
% image_mask, image_mask_norm
% label_fg, label_bg, label_mask
% pixel_mask, pixel_fg, pixel_bg
% idx_mask, idx_fg, idx_bg
% 盒子image_box，参数spacing，显示用的slice_show，存储用的label_save

% 在dataloader2的基础上将参数改成一个dict结构【2021.8】

image_path = dict_path.path_data_nii;
label_path_liver = dict_path.path_mask_nii;
% if(exist('dict_path.path_label_nii', 'var') && ~isempty(dict_path.path_label_nii)) % 判断结构里某个字段是否存在应该用isfield(dict_path, 'path_label_nii')
if(isfield(dict_path, 'path_label_nii') && ~isempty(dict_path.path_label_nii) && exist(dict_path.path_label_nii, 'file'))
    do_fg = true;
    label_path_vessel = dict_path.path_label_nii;
else
    do_fg = false;
end
idx_mask = dict_path.idx_mask;

%% 加载图像
image = load_nii(image_path);
spacing = image.hdr.dime.pixdim(2:4); % 参数spacing
image = image.img;
label_liver = load_nii(label_path_liver);
label_save = label_liver; % 留一个用来存
label_liver = label_liver.img;
if do_fg
    label_fg = load_nii(label_path_vessel);
    label_fg = label_fg.img;
else
    label_fg = [];
end

%% 【ATTENTION: 针对不同数据集要用不同的方式读取label】
% 贺建安的数据集中liver标注 = 6
% nnunet预分割结果liver标注 = 1
label_liver = label_liver==idx_mask;
% 我的数据集中有些vessel标注= 1，有些vessel标注肝静脉为1，门静脉为0
if do_fg
    label_fg = label_fg ~= 0;
    label_fg = logical(label_fg .* label_liver);
end

%% 计算box
image_box = GetBox3d(label_liver); % box要在去肿瘤之前算，不然不兼容

%% 裁剪
image_mask = image(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2));
% image_mask(label_mask) = 0;
label_mask = logical(label_liver(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2)));
if do_fg
    label_fg = logical(label_fg(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2)));
end

%% 图像加窗，40 350是腹腔脏器的窗，只是为了显示用！
WINDOW_POSITION = 40;
WINDOW_WIDTH = 350;
Imax = WINDOW_POSITION + WINDOW_WIDTH / 2;
Imin = WINDOW_POSITION - WINDOW_WIDTH / 2;
image_mask_norm = (double(image_mask) - Imin) / (Imax - Imin);

%% 用来显示的层数
slice_show = round((image_box(3,2) - image_box(3,1)) / 2);

end