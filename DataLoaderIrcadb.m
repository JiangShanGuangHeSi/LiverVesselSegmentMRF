function [image_mask, image_mask_norm, label_mask, image_box, spacing, slice_show, label_save, label_fg] ...
    = DataLoaderIrcadb(dict_path)
% 功能：读取ircadb nii格式的数据集，裁剪，返回
% image_mask, image_mask_norm
% label_fg, label_bg, label_mask
% pixel_mask, pixel_fg, pixel_bg
% idx_mask, idx_fg, idx_bg
% 盒子image_box，参数spacing，显示用的slice_show，存储用的label_save

% 将输入参数修改为dict结构，将返回值label_fg由第二个移动到最后一个【2021.8】

dataset_root = dict_path.root_dataset;
image_path = dict_path.path_image;
label_path_liver = dict_path.path_label_liver;
label_path_port = dict_path.path_label_port;
label_path_vena = dict_path.path_label_vena;


%% 加载图像
image = load_nii(image_path);
spacing = image.hdr.dime.pixdim(2:4); % 参数spacing
image = image.img;
label_liver = load_nii(label_path_liver);
label_save = label_liver; % 留一个用来存
label_liver = label_liver.img;
label_port = load_nii(label_path_port);
label_port = label_port.img;
label_vena = load_nii(label_path_vena);
label_vena = label_vena.img;

%% 计算box
image_box = GetBox3d(label_liver); % box要在去肿瘤之前算，不然不兼容

%% 补丁：去掉肿瘤。有些案例中（如case17）肿瘤非常影响结果
% 补丁：去掉case21中命名为Hyperplasie的肿瘤
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

%% 合并门静脉和肝静脉标签，限制在肝脏mask内
label_fg = label_liver ~= 0 & (label_port ~= 0 | label_vena ~= 0);% 肝脏内部血管标签索引
label_bg = label_liver ~= 0 & ~(label_port ~= 0 | label_vena ~= 0);% 肝脏内部背景（无血管）标签索引

%% 裁剪
label_fg = label_fg(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2));
label_bg = label_bg(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2));
label_mask = logical(label_liver(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2)));
image_mask = image(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2));
image_mask(~label_mask) = 0;

%% 图像加窗，40 350是腹腔脏器的窗，只是为了显示用！
WINDOW_POSITION = 40;
WINDOW_WIDTH = 350;
Imax = WINDOW_POSITION + WINDOW_WIDTH / 2;
Imin = WINDOW_POSITION - WINDOW_WIDTH / 2;
image_mask_norm = (image_mask - Imin) / (Imax - Imin);

%% 用来显示的层数
slice_show = round((image_box(3,2) - image_box(3,1)) / 2);

end