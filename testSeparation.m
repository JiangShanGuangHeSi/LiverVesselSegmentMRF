%% 测试大小血管分割
FOLDER_IMAGE = 'D:\\XueZhimeng\\project\\DataSet\\ZDYY\\ZDYYLiverVesselNII\\';
FOLDER_MASK = 'D:\\XueZhimeng\\project\\DataSet\\ZDYY\\ZDYYLiverPredictNnunet20210401\\';
FOLDER_LABEL = 'D:\\XueZhimeng\\project\\DataSet\\ZDYY\\ZDYYLiverVesselNIIlabel\\';
count = 0;
% record = struct(criterion,{}, parameters, {});
for idx = 1:59
    % 路径准备
    name_task = sprintf('ZDYY_%03d', idx);
    path_data_nii = [FOLDER_IMAGE, name_task, '_0000.nii.gz'];
    path_mask_nii = [FOLDER_MASK, name_task, '.nii.gz'];
    path_label_nii = [FOLDER_LABEL, name_task, '_label.nii.gz'];
    if ~exist(path_label_nii, 'file')
        continue
    else
        count = count + 1;
    end
    path_folder_process = ['result\\', name_task, '\\'];
    fprintf(['processing', name_task])
    
    [image_mask, image_mask_norm, label_mask, image_box, spacing, slice_show, label_save, label_fg, itpl]...
        = DataLoader4(path_data_nii, path_mask_nii, path_label_nii, 1);
    %%
    close all
    [label_large,label_small] = SeparationVesselSize(label_fg, 1);
    figure(1); 
    subplot(2, 3, 1); imshow(label_fg(:, :, slice_show));
    subplot(2, 3, 4); imshow(max(label_fg,[],3));
    subplot(2, 3, 2); imshow(label_large(:, :, slice_show));
    subplot(2, 3, 5); imshow(max(label_large,[],3));
    subplot(2, 3, 3); imshow(label_small(:, :, slice_show));
    subplot(2, 3, 6); imshow(max(label_small,[],3));
end

%% 测试不弹出窗口
handle=figure(1);
% set(handle, 'Visible', 'off');
set(handle, 'currentFigure', 0);