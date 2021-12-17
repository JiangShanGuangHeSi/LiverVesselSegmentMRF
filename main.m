% 批量处理中大医院和ircadb的数据
% 跑ircadb数据就注释掉 %% zdyy数据集 这一段代码
% 中大医院的数据注释掉 %% ircadb数据 这段代码

close all;
%% 将工具加入路径
addpath('NIfTI_20140122')
addpath('DPAD_toolbox')
addpath('Speckle-Reducing-Anisotropic-Diffusion-SRAD-master')
addpath('bm4d_matlab')
addpath('EnhancementFilter')
addpath('xzm_SymmetricConvexity\\SymmetricConvexityMEX1\\bin\\x64\\Debug')
addpath('phase-congruency-tensor-3d-master')

%% parameters: 参数
% bm4d参数：
parameters.bm4d.sigma             = 0.05;      % noise standard deviation given as percentage of the
                                               % maximum intensity of the signal, must be in [0,100]
parameters.bm4d.estimate_sigma    = 1;
parameters.bm4d.distribution_bm3d = 'Gauss'; % noise distribution
                                               %  'Gauss' --> Gaussian distribution
                                               %  'Rice ' --> Rician Distribution
parameters.bm4d.profile           = 'mp';    % BM4D parameter profile
                                               %  'lc' --> low complexity
                                               %  'np' --> normal profile
                                               %  'mp' --> modified profile
                                               % The modified profile is default in BM4D. For 
                                               % details refer to the 2013 TIP paper.
parameters.bm4d.do_wiener         = 1;       % Wiener filtering
                                               %  1 --> enable Wiener filtering
                                               %  0 --> disable Wiener filtering
parameters.bm4d.verbose           = 1;       % verbose mode
% sigmoid校正参数：
% parameters.sigmoid.correction_mode = 'sigmoid6';
parameters.sigmoid.correction_mode = 'sigmoid6';
parameters.sigmoid.kmeans_alpha = 15;
parameters.sigmoid.kmeans_beta = 125;
% 多尺度滤波参数：
parameters.enhancement_filter.scale = 1:6;
parameters.enhancement_filter.tau = 1;
% 相位一致性参数：
parameters.PC.sigmaOnf        = 0.01;%0.45    % Ratio of the standard deviation of the
                                              % Gaussian describing the log Gabor filter's
                                              % transfer function in the frequency domain
                                              % to the filter center frequency.        
parameters.PC.k               = 5;%5          % No of standard deviations of the noise
                                              % energy beyond the mean at which we set the
                                              % noise threshold point. 
parameters.PC.noiseMethod     = 0.2;%-1       % Choice of noise compensation method. 
% 多模型投票：
parameters.vote.mode = 'weight';
% parameters.vote.w_org = 0.6; 
% parameters.vote.w_vess = 0.3; 
% parameters.vote.w_PC = 0.1; 
% parameters.vote.thr_vote = 0.3;
parameters.vote.w_org = 0.7; 
parameters.vote.w_vess = 0.3; 
parameters.vote.w_PC = 0.0; 
parameters.vote.thr_vote = 0.4;
% ICM参数：
parameters.ICM.e = 2006;
parameters.ICM.Beta = [0.5 0.5];
parameters.ICM.IL = 1;
% parameters.ICM.weight_fg = 0.05;
parameters.ICM.weight_fg = -3;

%% ircadb数据
close all;
fid = fopen('data_path_info_ircadb.csv'); % 【ATTENTION: 这个文件里面的路径替换成你自己的】
Titles = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s',1,'delimiter', ',');
Pathes = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %d','delimiter', ',');
fclose(fid);

dict_path = {};
for path_idx = 1:22
    close all;
    name_task = sprintf('ircadb_%03d', path_idx);
    path_folder_process = ['result/', name_task, '/'];
    path_folder_featuremap = ['feature_map\\', name_task, '\\'];
    %% 加载/裁剪图像
    dict_path.root_dataset = ['D:\XueZhimeng\project\DataSet\3Dircadb\NiiFile\case', char(num2str(path_idx)), '\']; % 【ATTENTION: 这个路径替换成你自己的】
    dict_path.path_image = Pathes{1,1}{path_idx,1};
    dict_path.path_label_liver = Pathes{1,2}{path_idx,1};
    dict_path.path_label_port = Pathes{1,3}{path_idx,1};
    dict_path.path_label_vena = Pathes{1,4}{path_idx,1};
    if path_idx == 1
        record = MLE_MRF5_forAllDataset(@DataLoaderIrcadb, dict_path, path_folder_process, parameters, false, path_folder_featuremap);
    else
        record(path_idx) = MLE_MRF5_forAllDataset(@DataLoaderIrcadb, dict_path, path_folder_process, parameters, false, path_folder_featuremap);
    end
end
count = 22;

%% zdyy数据集
% FOLDER_IMAGE = 'D:\\XueZhimeng\\project\\DataSet\\ZDYY\\ZDYYLiverVesselNII\\'; % 【ATTENTION: 这个路径替换成你自己的】
% FOLDER_MASK = 'D:\\XueZhimeng\\project\\DataSet\\ZDYY\\ZDYYLiverPredictNnunet20210401\\'; % 【ATTENTION: 这个路径替换成你自己的】
% FOLDER_LABEL = 'D:\\XueZhimeng\\project\\DataSet\\ZDYY\\ZDYYLiverVesselNIIlabel\\'; % 【ATTENTION: 这个路径替换成你自己的】
% count = 0;
% % record = struct(criterion,{}, parameters, {});
% dict_path = {};
% record = {};
% zdyy_idx = [];
% for idx = 1:59 %[16, 49, 58]%
%     % 路径准备
%     name_task = sprintf('ZDYY_%03d', idx);
%     dict_path.path_data_nii = [FOLDER_IMAGE, name_task, '_0000.nii.gz'];
%     dict_path.path_mask_nii = [FOLDER_MASK, name_task, '.nii.gz'];
%     dict_path.path_label_nii = [FOLDER_LABEL, name_task, '_label.nii.gz'];
%     dict_path.idx_mask = 1;
%     if ~exist(dict_path.path_label_nii, 'file')
%         continue
%     else
%         count = count + 1;
%         zdyy_idx = [zdyy_idx; idx];
%     end
%     path_folder_process = ['result\\', name_task, '\\'];
%     fprintf(['processing', name_task])
%     
%     path_folder_featuremap = ['feature_map\\', name_task, '\\'];
%     
%     % 计算
% %     MLE_MRF1(path_data_nii, path_mask_nii, path_folder_process, parameters);
%     if count == 1
% %         record = MLE_MRF2(path_data_nii, path_mask_nii, path_label_nii, path_folder_process, parameters, false, path_folder_featuremap);
%         record = MLE_MRF5_forAllDataset(@DataLoaderZDYY, dict_path, path_folder_process, parameters, false, path_folder_featuremap);
%     else
% %         record(count) = MLE_MRF2(path_data_nii, path_mask_nii, path_label_nii, path_folder_process, parameters, false, path_folder_featuremap);
%         record(count) = MLE_MRF5_forAllDataset(@DataLoaderZDYY, dict_path, path_folder_process, parameters, false, path_folder_featuremap);
%     end
% %     if count == 1
% %         record = MLE_MRF2_interpPC(path_data_nii, path_mask_nii, path_label_nii, path_folder_process, parameters);
% %     else
% %         record(count) = MLE_MRF2_interpPC(path_data_nii, path_mask_nii, path_label_nii, path_folder_process, parameters);
% %     end
% end

%% 解析record结果
% 编号·不同特征图·不同指标 表
n_Featuremap = length(fieldnames(record(1).criterion.Dx_vote));
n_Criterion = length(fieldnames(record(1).criterion));
% crit = zeros(count, n_Featuremap, n_Criterion);
crit = zeros(count, n_Criterion, n_Featuremap);
for idx = 1:count
    crit(idx, 1, :) = cell2mat(struct2cell(record(idx).criterion.MLE_vess));
    crit(idx, 2, :) = cell2mat(struct2cell(record(idx).criterion.MLE_PC));
    crit(idx, 3, :) = cell2mat(struct2cell(record(idx).criterion.MLE_org));
    crit(idx, 4, :) = cell2mat(struct2cell(record(idx).criterion.Dx_vote));
    crit(idx, 5, :) = cell2mat(struct2cell(record(idx).criterion.Dx_vote_removeisland));
    crit(idx, 6, :) = cell2mat(struct2cell(record(idx).criterion.ICM1));
    % 4
%     crit(idx, 7, :) = cell2mat(struct2cell(record(idx).criterion.ICM2));
%     crit(idx, 8, :) = cell2mat(struct2cell(record(idx).criterion.ICM3));
%     crit(idx, 9, :) = cell2mat(struct2cell(record(idx).criterion.ICM4));
%     crit(idx, 10, :) = cell2mat(struct2cell(record(idx).criterion.ICM_removeisland));
%     crit(idx, 11, :) = cell2mat(struct2cell(record(idx).criterion.ICM_large));
%     crit(idx, 12, :) = cell2mat(struct2cell(record(idx).criterion.ICM_small));
%     crit(idx, 13, :) = cell2mat(struct2cell(record(idx).criterion.ICM_large2orglabel));
    % 1
    crit(idx, 7, :) = cell2mat(struct2cell(record(idx).criterion.ICM_removeisland));
    crit(idx, 8, :) = cell2mat(struct2cell(record(idx).criterion.ICM_large));
    crit(idx, 9, :) = cell2mat(struct2cell(record(idx).criterion.ICM_small));
    crit(idx, 10, :) = cell2mat(struct2cell(record(idx).criterion.ICM_large2orglabel));
%     % 2
%     crit(idx, 7, :) = cell2mat(struct2cell(record(idx).criterion.ICM2));
%     crit(idx, 8, :) = cell2mat(struct2cell(record(idx).criterion.ICM_removeisland));
%     crit(idx, 9, :) = cell2mat(struct2cell(record(idx).criterion.ICM_large));
%     crit(idx, 10, :) = cell2mat(struct2cell(record(idx).criterion.ICM_small));
%     crit(idx, 11, :) = cell2mat(struct2cell(record(idx).criterion.ICM_large2orglabel));
end
crit_acc = crit(:,:,1);
crit_sens = crit(:,:,2);
crit_spec = crit(:,:,3);
crit_prec = crit(:,:,4);
crit_dice = crit(:,:,5);
crit_loU = crit(:,:,6);
crit1 = crit;
crit1([21, 22], :, :) = [];
% crit1([3, 9, 21, 22], :, :) = [];
% crit1([3, 21, 22], :, :) = [];
% 去掉7,42,58
% crit1([3, 15, 23], :, :) = [];
% crit1([3, 15, 20, 22, 23], :, :) = [];
crit_all = squeeze(mean(crit1, 1));