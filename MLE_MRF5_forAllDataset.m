function [record] = MLE_MRF5_forAllDataset(dataloader, dict_path_data, path_folder_process, parameters, fig_visible, path_folder_featuremap)
% ������MLE_MRF1�������˿�ѡ��label��
% ��MLE��MRF��������ͼ��ȫ����
% path_data_nii: ԭʼ����nii�ļ�·��
% path_mask_nii: ������ģnii�ļ�·��
% ��������path_label_nii: ��ѡ�������label��������ִ�����ۣ����Ϊ�����Ԥ�⡿
% ����������ֵ��record������ʵ���ϵĲ���parameters�����۽��criterion��
% ����������
% ��������ԭʼͼ��ȴ������ͼ.m�ı��棬�����˸���ͼ.nii�ļ����桿
% ��Ϊ�˼���ircadb�������ݼ���·���ϲ�Ϊһ����dict_path_data��������Ҫ���ݶ�Ӧdataloader�ĺ��������
% ircadb���ݼ���dict_path_dataӦ�ð�����dataset_root, image_path, label_path_liver, label_path_port, label_path_vena
% �д�ҽԺ���ݼ���dict_path_dataӦ�ð�����path_data_nii, path_mask_nii, path_label_nii(��ѡ), idx_mask

% path_folder_process: �м���̱���·��������
%   ������parameters.mat
%   ����ͼ��image_bm4d.mat, image_vesselness.mat, image_CP.mat, 
%   ����ͼ��image_Dp_org.mat, image_Dp_vess.mat, image_Dp_CP.mat
%   �����parameters.IT�� ICMx.nii.gzͼ��
%   �Լ������м����ͼ��

% % parameters: ����
% % bm4d������
%     parameters.bm4d.sigma             = 0.05;      % noise standard deviation given as percentage of the
%                                                    % maximum intensity of the signal, must be in [0,100]
%     parameters.bm4d.estimate_sigma    = 1;
%     parameters.bm4d.distribution_bm3d = 'Gauss'; % noise distribution
%                                                    %  'Gauss' --> Gaussian distribution
%                                                    %  'Rice ' --> Rician Distribution
%     parameters.bm4d.profile           = 'mp';    % BM4D parameter profile
%                                                    %  'lc' --> low complexity
%                                                    %  'np' --> normal profile
%                                                    %  'mp' --> modified profile
%                                                    % The modified profile is default in BM4D. For 
%                                                    % details refer to the 2013 TIP paper.
%     parameters.bm4d.do_wiener         = 1;       % Wiener filtering
%                                                    %  1 --> enable Wiener filtering
%                                                    %  0 --> disable Wiener filtering
%     parameters.bm4d.verbose           = 1;       % verbose mode
% % sigmoidУ��������
%     parameters.sigmoid.correction_mode
%     parameters.sigmoid.kmeans_alpha = 15;
%     parameters.sigmoid.kmeans_beta = 125;
% % ��߶��˲�������
%     parameters.enhancement_filter.scale = 1:6;
%     parameters.enhancement_filter.tau = 1;
% % ��λһ���Բ�����
%     parameters.PC.sigmaOnf        = 0.01;%0.45    % Ratio of the standard deviation of the
%                                                   % Gaussian describing the log Gabor filter's
%                                                   % transfer function in the frequency domain
%                                                   % to the filter center frequency.        
%     parameters.PC.k               = 5;%5          % No of standard deviations of the noise
%                                                   % energy beyond the mean at which we set the
%                                                   % noise threshold point. 
%     parameters.PC.noiseMethod     = 0.2;%-1       % Choice of noise compensation method. 
% % ��ģ��ͶƱ��
%     parameters.vote.mode = 'weight';
%     parameters.vote.w_org = 0.6; 
%     parameters.vote.w_vess = 0.3; 
%     parameters.vote.w_PC = 0.1; 
%     parameters.vote.thr_vote = 0.3;
% % ICM������
%     parameters.ICM.e = 2006;
%     parameters.ICM.Beta = [1 1];
%     parameters.ICM.IL = 4;
%     parameters.ICM.weight_fg = 0.05;
if nargin < 6
    fig_visible = false; %Ĭ�ϲ���ʾͼ��
end

close all;

% ͼ��Ӵ���40 350�Ǹ�ǻ�����Ĵ���ֻ��Ϊ����ʾ�ã�
WINDOW_POSITION = 40;
WINDOW_WIDTH = 350;
Imax = WINDOW_POSITION + WINDOW_WIDTH / 2;
Imin = WINDOW_POSITION - WINDOW_WIDTH / 2;

%% �����û�оͣ�����һ�����ļ���
[status,msg,msgID] = mkdir(path_folder_process);
[status,msg,msgID] = mkdir(path_folder_featuremap);

%% ����ͼ���ֲ���ȷ�ϣ���һ����Ϊ��ȡ����ȥ�Ľ�����Ӷ���Լʱ�䣩
% % �ж��Ƿ����label����MLE_MRF5�汾�����ˣ�Ĭ�϶����ڣ���Ȼdataloaderû��Ū������
% if exist('path_label_nii', 'var') && ~isempty(path_label_nii)
%     do_crit = true;
% else
%     do_crit = false;
% end
% do_crit = true;
% ������оͣ������ϴβ�����ȷ����Щ�����Ƿ�Ҫ����
if exist([path_folder_process, 'parameters.mat'])
    parameters_old = load([path_folder_process, 'parameters.mat']);
    parameters_old = parameters_old.parameters;
    % ȷ���Ƿ����BM4D
    if (parameters_old.bm4d.sigma == parameters.bm4d.sigma && ...
            parameters_old.bm4d.estimate_sigma == parameters.bm4d.estimate_sigma && ...
            strcmp(parameters_old.bm4d.distribution_bm3d, parameters.bm4d.distribution_bm3d) && ...
            strcmp(parameters_old.bm4d.profile, parameters.bm4d.profile) && ...
            parameters_old.bm4d.do_wiener == parameters.bm4d.do_wiener && ...
            parameters_old.bm4d.verbose == parameters.bm4d.verbose)
        do_bm4d = false;
    else
        do_bm4d = true;
    end
    % ȷ���Ƿ���ж�߶��˲�
    if (strcmp(parameters_old.sigmoid.correction_mode, parameters.sigmoid.correction_mode) &&...
            all(parameters_old.enhancement_filter.scale == parameters.enhancement_filter.scale) && ...
            parameters_old.enhancement_filter.tau == parameters.enhancement_filter.tau)
        do_enhancement_filter = false;
        if(strcmp(parameters.sigmoid.correction_mode, 'mannual sigmoid') && ...
                (parameters.sigmoid.kmeans_alpha ~= parameters.sigmoid.kmeans_alpha || ...
                parameters.sigmoid.kmeans_beta ~= parameters.sigmoid.kmeans_beta))
            do_enhancement_filter = true;
        end
    else
        do_enhancement_filter = true;
    end
    % ȷ���Ƿ������λһ����
    if (parameters_old.PC.sigmaOnf == parameters.PC.sigmaOnf && ...
            parameters_old.PC.k == parameters.PC.k && ...
            parameters_old.PC.noiseMethod == parameters.PC.noiseMethod)
        do_PC = false;
    else
        do_PC = true;
    end
else % ���û��parameters����ʾ֮ǰû�н��й����㣬�������
    do_bm4d = true;
    do_enhancement_filter = true;
    do_PC = true;
end
% do_enhancement_filter = true;
% do_PC = true; % ����
%% �������
save( [path_folder_process, 'parameters.mat'], 'parameters');

%% ����ͼ��
[image_mask, image_norm, label_mask, image_box, spacing, slice_show, label_save, label_fg] ...
    = dataloader(dict_path_data);
if exist('label_fg', 'var') && ~isempty(label_fg)
    do_crit = true;
else
    do_crit = false;
end
% ������label����weight��
% parameters.ICM.weight_fg = sum(label_fg(:)) / sum(label_mask(:));
% ������ڣ������ǩͶӰ
if do_crit
imwrite(label_fg(:, :, slice_show), [path_folder_process, 'image_0_label.png'])
imwrite(max(label_fg,[],3), [path_folder_process, 'image_0_label_projection.png'])
end

%% BM4D
dir_image_bm4d = [path_folder_process, 'image_bm4d.mat'];
if exist(dir_image_bm4d, 'file') && ~do_bm4d
    load(dir_image_bm4d);
else
    [image_bm4d, image_bm4d_sigma_est] = bm4d(image_norm, ...
        parameters.bm4d.distribution_bm3d, ...
        (~parameters.bm4d.estimate_sigma)*parameters.bm4d.sigma, ...
        parameters.bm4d.profile, parameters.bm4d.do_wiener, parameters.bm4d.verbose);
    save(dir_image_bm4d, 'image_bm4d', 'image_bm4d_sigma_est');
end
image_bm4d = image_bm4d.*(label_mask~=0);
imwrite(image_norm(:, :, slice_show), [path_folder_process, 'image.png'])
imwrite(image_bm4d(:, :, slice_show), [path_folder_process, 'image_1_bm4d.png'])
% ����ԭʼͼ���С��ͼƬ.mat
image_save = zeros(size(label_save.img));
image_save(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2))...
        = image_bm4d;
save( [path_folder_featuremap, 'image_bm4d_s.mat'], 'image_save');

%% EM�㷨׼��
image_EM = round((image_bm4d*(Imax-Imin)+Imin));
pixels_EM = image_EM(label_mask);    % bm4d
hist_EM = GetDecomposedHist(round(image_EM), [], label_mask);
%     pixels_EM = pixels_liver_bm4d_mip;    % bm4d mip
%     hist_EM = hist_bm4d_mip;
%     image_EM = round((image_bm4d_mip*(Imax-Imin)+Imin));

range = min(pixels_EM):max(pixels_EM);
image_EM_norm = (image_EM - min(image_EM(:))) / (max(image_EM(:) - min(image_EM(:))));

%% sigmoidУ��
if parameters.sigmoid.correction_mode == "sigmoid6"
    % kmeans�����ҳ�sigmoid grayscale mapping����alpha_1��beta_1
    % ��ʼ������ֱ�Ϊǰ0.05, 0.20, 0.5, 0.7, 0.95��λ��
%     [~, position] = min(abs(cumsum(hist_EM.N_mask') / sum(hist_EM.N_mask) - [0.08, 0.3, 0.5, 0.7, 0.92]), [], 1);
    [~, position] = min(abs(cumsum(hist_EM.N_mask') / sum(hist_EM.N_mask) - [0.2, 0.3, 0.7, 0.93, 0.95, 0.97]), [], 1);
%     [kmeans_idx, kmeans_Ctrs] = kmeans(pixels_EM, 5, 'start', transpose(hist_EM.X_mask(position))); % The 'Start' matrix must have K rows.
    [kmeans_idx, kmeans_Ctrs] = kmeans(pixels_EM, 6, 'start', transpose(hist_EM.X_mask(position))); % The 'Start' matrix must have K rows.
    [kmeans_idx, kmeans_Ctrs] = Kmean_reorder(kmeans_idx,kmeans_Ctrs);

    % sigmoid У��
%     kmeans_alpha = (kmeans_Ctrs(5) - kmeans_Ctrs(4)) / 2; % ��5��
%     kmeans_beta = (kmeans_Ctrs(5) + kmeans_Ctrs(4)) / 2;
    parameters.sigmoid.kmeans_alpha = (kmeans_Ctrs(6) - kmeans_Ctrs(5)) / 2;% ��6��
%     kmeans_beta = (kmeans_Ctrs(6) + kmeans_Ctrs(5)) / 2;
    parameters.sigmoid.kmeans_beta = (kmeans_Ctrs(6) + kmeans_Ctrs(5)) / 2 + 10; % ��У��
    image_sigmoid = 1 / (1 + exp((parameters.sigmoid.kmeans_beta - image_EM) / parameters.sigmoid.kmeans_alpha)); % ����
elseif parameters.sigmoid.correction_mode == "sigmoid5"
    % kmeans�����ҳ�sigmoid grayscale mapping����alpha_1��beta_1
    % ��ʼ������ֱ�Ϊǰ0.05, 0.20, 0.5, 0.7, 0.95��λ��
%     [~, position] = min(abs(cumsum(hist_EM.N_mask') / sum(hist_EM.N_mask) - [0.08, 0.3, 0.5, 0.7, 0.92]), [], 1);
    [~, position] = min(abs(cumsum(hist_EM.N_mask') / sum(hist_EM.N_mask) - [0.2, 0.3, 0.7, 0.93, 0.95]), [], 1);
%     [kmeans_idx, kmeans_Ctrs] = kmeans(pixels_EM, 5, 'start', transpose(hist_EM.X_mask(position))); % The 'Start' matrix must have K rows.
    [kmeans_idx, kmeans_Ctrs] = kmeans(pixels_EM, 5, 'start', transpose(hist_EM.X_mask(position))); % The 'Start' matrix must have K rows.
    [kmeans_idx, kmeans_Ctrs] = Kmean_reorder(kmeans_idx,kmeans_Ctrs);

    % sigmoid У��
    parameters.sigmoid.kmeans_alpha = (kmeans_Ctrs(5) - kmeans_Ctrs(4)) / 2; % ��5��
    parameters.sigmoid.kmeans_beta = (kmeans_Ctrs(5) + kmeans_Ctrs(4)) / 2;
    image_sigmoid = 1 / (1 + exp((parameters.sigmoid.kmeans_beta - image_EM) / parameters.sigmoid.kmeans_alpha)); % ����
elseif parameters.sigmoid.correction_mode == "mannual linear"
%         WINDOW_WIDTH = [80, 50, 90, 115, 110, 100, 70, 100, 80, 80,...
%             70, 50, 70, 90, 90, 100, 80, 90, 70, 80, 70];
%         WINDOW_LEVEL = [170, 140, 130, 185, 160, 190, 160, 250, 150, 180,...
%             160, 170, 190, 170, 180, 120, 160, 120, 140, 110, 140];
%������ͼʡ��û��д����
    WINDOW_WIDTH = [70, 90, 150, 80, 100, 90, 70, 100, 80, 80,...
        70, 50, 70, 90, 90, 100, 80, 90, 70, 80, 70];
    WINDOW_LEVEL = [170, 150, 180, 140, 160, 160, 160, 250, 150, 180,...
        160, 170, 190, 170, 180, 120, 160, 120, 140, 110, 140];
    Imax_correction = WINDOW_LEVEL(path_idx) + WINDOW_WIDTH(path_idx) / 2;
    Imin_correction = WINDOW_LEVEL(path_idx) - WINDOW_WIDTH(path_idx) / 2;
%         image_sigmoid = (image_liver - Imin_correction) / (Imax_correction - Imin_correction);
    image_sigmoid = (image_EM - Imin_correction) / (Imax_correction - Imin_correction);
elseif parameters.sigmoid.correction_mode == "mannual sigmoid"
    image_sigmoid = 1 / (1 + exp((parameters.sigmoid.kmeans_beta - image_EM) / parameters.sigmoid.kmeans_alpha)); % ����
end

figure(5); fig5 = gcf;
if fig_visible == false
    set(fig5, 'Visible', 'off'); % ������
end
set(fig5,'Units','centimeter','Position',[5 5 12 15]);% ���ô��ڴ�С
subplot(3,3,1); imshow((image_EM(:, :, slice_show)-Imin)/(Imax-Imin)); title('image denoise');
% subplot(3,3,2); imshow(label_fg(:,:, slice_show)); title('label');
% subplot(3,3,3); imshow(max(label_fg,[],3)); title('label projection');
subplot(3,3,4); imshow(image_sigmoid(:, :, slice_show)); title('sigmoid');

imwrite(image_sigmoid(:, :, slice_show), [path_folder_process, 'image_1_sigmoid.png'])
imwrite(max(image_sigmoid,[],3), [path_folder_process, 'image_1_sigmoid_projection.png'])

% ����ԭʼͼ���С��ͼƬ.mat
image_save = zeros(size(label_save.img));
image_save(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2))...
        = image_sigmoid;
save( [path_folder_featuremap, 'image_sigmoid_s.mat'], 'image_save');

%% ��߶��˲���spacingѡ����ʵ��λmm�����scaleҲҪѡ��Ѫ�ܵ�λΪmm�˳߶ȷ�Χ
dir_image_vesselness = [path_folder_process, 'image_vesselness.mat'];
if exist(dir_image_vesselness, 'file') && ~do_enhancement_filter
    load(dir_image_vesselness);
else
    [image_vesselness, Vx, Vy, Vz] = vesselness3D(image_sigmoid, parameters.enhancement_filter.scale, spacing, parameters.enhancement_filter.tau, true);
    % Zhang�ķ���
%     [image_vesselness, Vx, Vy, Vz] = vesselness3DwithVectorZhang(image_sigmoid, parameters.enhancement_filter.scale, spacing, parameters.enhancement_filter.tau, true);
    save(dir_image_vesselness, 'image_vesselness', 'Vx', 'Vy', 'Vz');
end
subplot(3,3,5); imshow(image_vesselness(:,:, slice_show)); title('vesselness');
subplot(3,3,6); imshow(max(image_vesselness,[],3)); title('vesselness projection');

imwrite(image_vesselness(:,:, slice_show), [path_folder_process, 'image_1_vess.png']);
imwrite(max(image_vesselness,[],3), [path_folder_process, 'image_1_vess_projection.png']);

% ����ԭʼͼ���С��ͼƬ.mat
image_save = zeros(size(label_save.img));
image_save(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2))...
        = image_vesselness;
save( [path_folder_featuremap, 'image_vesselness_s.mat'], 'image_save');

%% ��λһ����
dir_image_PC = [path_folder_process, 'image_PC.mat'];
if exist(dir_image_PC, 'file') && ~do_PC
    load(dir_image_PC);
else
%     [image_PC,image_PV,image_PN] = GetPhaseCongruency3D(image_sigmoid, parameters.PC.sigmaOnf, ...
%         parameters.PC.k, parameters.PC.noiseMethod);% �Ƿ�ת�ƺ��Խ��û��Ӱ��
    [image_PC,image_PV,image_PN] = GetPhaseCongruency3D(image_sigmoid);% �Ƿ�ת�ƺ��Խ��û��Ӱ��
    save(dir_image_PC, 'image_PC', 'image_PV', 'image_PN');
end
figure(6); fig6 = gcf;
if fig_visible == false
    set(fig6, 'Visible', 'off'); % ������
end
set(fig6,'Units','centimeter','Position',[5 5 12 12]);% ���ô��ڴ�С
subplot(2,2,1); imshow(max(image_sigmoid,[],3)); title('sigmoid projection');
subplot(2,2,2); imshow(max(image_PC,[],3)); title('Phase Congruency projection');
subplot(2,2,3); imshow(max(image_PV / max(image_PV(:)),[],3)); title('Phase Vesselness projection');
subplot(2,2,4); imshow(max(image_PN,[],3)); title('Phase Neuritenees projection');
figure(8); fig8 = gcf;
if fig_visible == false
    set(fig8, 'Visible', 'off'); % ������
end
set(fig8,'Units','centimeter','Position',[5 5 12 12]);% ���ô��ڴ�С
subplot(2,2,1); imshow(image_sigmoid(:, :, slice_show)); title('sigmoid projection');
subplot(2,2,2); imshow(image_PC(:, :, slice_show)); title('Phase Congruency projection');
subplot(2,2,3); imshow(image_PV(:, :, slice_show) / max(image_PV(:))); title('Phase Vesselness projection');
subplot(2,2,4); imshow(image_PN(:, :, slice_show)); title('Phase Neuritenees projection');

imwrite(image_PC(:,:, slice_show), [path_folder_process, 'image_1_PC.png']);
imwrite(max(image_PC,[],3), [path_folder_process, 'image_1_PC_projection.png']);

% ����ԭʼͼ���С��ͼƬ.mat
image_save = zeros(size(label_save.img));
image_save(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2))...
        = image_PC;
save( [path_folder_featuremap, 'image_PC_s.mat'], 'image_save');

%% ���ʷֲ�ģ�ͣ���ˮƽģ�ͣ�
[MLE_vess.distribution, MLE_vess.Object, MLE_vess.VBN_Init, MLE_vess.VBN_EM, MLE_vess.Dx_MLE, MLE_vess.sort_D, MLE_vess.image_Dp] ...
    = DistributionInitialize_EM_MLE('eeg', [0.1; 0.3; 0.6], image_vesselness, label_mask, [], 1, 1);
[MLE_PC.distribution, MLE_PC.Object, MLE_PC.VBN_Init, MLE_PC.VBN_EM, MLE_PC.Dx_MLE, MLE_PC.sort_D, MLE_PC.image_Dp] ...
    = DistributionInitialize_EM_MLE('eg', [0.1; 0.6], image_PC, label_mask, [], 1, 1);
[MLE_org.distribution, MLE_org.Object, MLE_org.VBN_Init, MLE_org.VBN_EM, MLE_org.Dx_MLE, MLE_org.sort_D, MLE_org.image_Dp] ...
    = DistributionInitialize_EM_MLE('ggg', [], round((image_bm4d*(Imax-Imin)+Imin)), label_mask, MLE_vess.Dx_MLE, 10, 3);

imwrite(MLE_vess.image_Dp(:,:, slice_show), [path_folder_process, 'image_2_vess_Dp.png']);
imwrite(max(MLE_vess.image_Dp,[],3), [path_folder_process, 'image_2_vess_Dp_projection.png']);
imwrite(MLE_PC.image_Dp(:,:, slice_show), [path_folder_process, 'image_2_PC_Dp.png']);
imwrite(max(MLE_PC.image_Dp,[],3), [path_folder_process, 'image_2_PC_Dp_projection.png']);
imwrite(MLE_org.image_Dp(:,:, slice_show), [path_folder_process, 'image_2_org_Dp.png']);
imwrite(max(MLE_org.image_Dp,[],3), [path_folder_process, 'image_2_org_Dp_projection.png']);

% ����ԭʼͼ���С��ͼƬ.mat
image_save = zeros(size(label_save.img));
image_save(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2))...
        = MLE_vess.image_Dp;
save( [path_folder_featuremap, 'image_vesselness_MLE_s.mat'], 'image_save');
image_save = zeros(size(label_save.img));
image_save(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2))...
        = MLE_PC.image_Dp;
save( [path_folder_featuremap, 'image_PC_MLE_s.mat'], 'image_save');
image_save = zeros(size(label_save.img));
image_save(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2))...
        = MLE_org.image_Dp;
save( [path_folder_featuremap, 'image_org_MLE_s.mat'], 'image_save');

% ����ԭʼͼ���С��nii
label_predict = zeros(size(label_save.img), 'uint8');
label_predict(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2))...
    = uint8(MLE_vess.image_Dp * 1000);
label_save.img = label_predict;
save_nii(label_save, [path_folder_featuremap, 'image_vesselness_MLE_s.nii.gz']);
label_predict = zeros(size(label_save.img), 'uint8');
label_predict(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2))...
    = uint8(MLE_PC.image_Dp * 1000);
label_save.img = label_predict;
save_nii(label_save, [path_folder_featuremap, 'image_PC_MLE_s.nii.gz']);
label_predict = zeros(size(label_save.img), 'uint8');
label_predict(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2))...
    = uint8(MLE_org.image_Dp * 1000);
label_save.img = label_predict;
save_nii(label_save, [path_folder_featuremap, 'image_org_MLE_s.nii.gz']);

%% ����record
record.parameters = parameters;
% record.criterion = crit;
% return %����

%% ���۸��ʷֲ�ģ�͵�Ч��
if do_crit
    fprintf('criterion MLE vess: ')
    [crit.MLE_vess.accuracy, crit.MLE_vess.sensitivity, crit.MLE_vess.specificity, crit.MLE_vess.precision, ...
        crit.MLE_vess.dice, crit.MLE_vess.loU] = criterion(MLE_vess.Dx_MLE, label_fg);
    fprintf('criterion MLE PC: ')
    [crit.MLE_PC.accuracy, crit.MLE_PC.sensitivity, crit.MLE_PC.specificity, crit.MLE_PC.precision, ...
        crit.MLE_PC.dice, crit.MLE_PC.loU] = criterion(MLE_PC.Dx_MLE, label_fg);
    fprintf('criterion MLE org: ')
    [crit.MLE_org.accuracy, crit.MLE_org.sensitivity, crit.MLE_org.specificity, crit.MLE_org.precision, ...
        crit.MLE_org.dice, crit.MLE_org.loU] = criterion(MLE_org.Dx_MLE, label_fg);
end

%% ͶƱ
if parameters.vote.mode == "majority" % ����ͶƱ
    Dx_vote = Object * (((MLE_org.Dx_MLE>0) + (MLE_vess.Dx_MLE>0) + (MLE_PC.Dx_MLE>0)) >= 2);
elseif parameters.vote.mode == "weight" % ��ȨͶƱ
    image_Dp_vote = parameters.vote.w_org * MLE_org.image_Dp ...
        + parameters.vote.w_vess * MLE_vess.image_Dp + ...
        parameters.vote.w_PC * MLE_PC.image_Dp;
    Dx_vote = image_Dp_vote > parameters.vote.thr_vote;
    
    imwrite(image_Dp_vote(:,:, slice_show), [path_folder_process, 'image_3_vote_Dp.png']);
    imwrite(max(image_Dp_vote,[],3), [path_folder_process, 'image_3_vote_Dp_projection.png']);
end
% Dx_vote_removeisland = MLE_org.Object * RemoveIsland(Dx_vote, 50);
Dx_vote_removeisland = MLE_org.Object * RemoveIsland(Dx_vote, 100); %��20210817��
figure(9); fig9 = gcf;
if fig_visible == false
    set(fig9, 'Visible', 'off'); % ������
end
set(fig9,'Units','centimeter','Position',[5 5 12 12]);% ���ô��ڴ�С
subplot(2,2,1); imshow(max(MLE_org.image_Dp,[],3)); title('Dp org projection');
subplot(2,2,2); imshow(max(MLE_vess.image_Dp,[],3)); title('Dp vess projection');
subplot(2,2,3); imshow(max(MLE_PC.image_Dp,[],3)); title('Dp PC projection');
subplot(2,2,4); imshow(max(image_Dp_vote,[],3)); title('Dp vote projection');

imwrite(Dx_vote(:,:, slice_show) * 255, [path_folder_process, 'image_3_vote_Dx.png']);
imwrite(max(Dx_vote,[],3) * 255, [path_folder_process, 'image_3_vote_Dx_projection.png']);

%% ����ͶƱ���
if do_crit
    fprintf('criterion MLE vote: ')
    [crit.Dx_vote.accuracy, crit.Dx_vote.sensitivity, crit.Dx_vote.specificity, crit.Dx_vote.precision, ...
        crit.Dx_vote.dice, crit.Dx_vote.loU] = criterion(Dx_vote, label_fg);
    fprintf('criterion MLE vote remove island: ')
    [crit.Dx_vote_removeisland.accuracy, crit.Dx_vote_removeisland.sensitivity, crit.Dx_vote_removeisland.specificity, crit.Dx_vote_removeisland.precision, ...
        crit.Dx_vote_removeisland.dice, crit.Dx_vote_removeisland.loU] = criterion(Dx_vote_removeisland, label_fg);
end


%% ICM
% parameters.ICM.e=2013;% ���e����Ӧ����3�������Ǹ�������ûʲô����
% Dx = ICM_MRF(VBN_EM,flagIMG,Dx_MLE,sort_D,image_liver,image_vesselness,Object,NB,IL,Beta,Vx,Vy,Vz,e,slice_show);
inPLX = MLE_org.sort_D(:,MLE_org.Object+1);
pxl_k = MLE_org.sort_D(:,1:MLE_org.Object);
% inPLX = MLE_PC.sort_D(:,MLE_PC.Object+1); % ������Ϲ�㡿
% pxl_k = MLE_PC.sort_D(:,1:MLE_PC.Object);
Dx_init = Dx_vote_removeisland; % Dx_MLE_removeisland;
OptiumBeTa = parameters.ICM.Beta; % �ⲻ�Ż�beta����˵û��
figure(14)
fig14 = gcf;
if fig_visible == false
    set(fig14, 'Visible', 'off'); % ������
end
set(fig14,'Units','centimeter','Position',[5 5 20 5]);
D0 = zeros(size(Dx_init)); % ׼��D0�ռ�

if parameters.ICM.weight_fg > 0
    weight_fg = parameters.ICM.weight_fg;
elseif parameters.ICM.weight_fg == -1
    weight_fg = sum(MLE_org.Dx_MLE(:)) / sum(label_mask(:));
elseif parameters.ICM.weight_fg == -2
    weight_fg = sum(Dx_vote(:) > 0) / sum(label_mask(:));
elseif parameters.ICM.weight_fg == -3
    weight_fg = sum(label_fg(:)) / sum(label_mask(:)); % ���ף�
end

% weight_fg = sum(Dx_vote(:) > 0) / sum(label_mask(:));
% weight_fg = sum(label_fg(:)) / sum(label_mask(:));
for t = 1:parameters.ICM.IL
%     tic;
    disp(['BeTa = ' num2str(OptiumBeTa) ';' 'ICM iteration ' num2str(t) ' times...']); 
%     weight_fg = sum(label_fg(:)) / sum(label_mask(:));
    [Dx, image_pV, image_pB] = ICM_estimation(MLE_org.VBN_EM,pxl_k,inPLX,MLE_org.Object,Dx_init,OptiumBeTa,Vx,Vy,Vz,image_EM,parameters.ICM.e,weight_fg);
%     % ��Ϲ��MLE_PC.Object��
%     [Dx, image_pV, image_pB] = ICM_estimation(MLE_org.VBN_EM,pxl_k,inPLX,MLE_PC.Object,Dx_init,OptiumBeTa,Vx,Vy,Vz,image_EM,parameters.ICM.e,weight_fg);

    figure(14)
    if fig_visible == false
        set(fig14, 'Visible', 'off'); % ������
    end
    subplot(max(1,floor(parameters.ICM.IL/5)),5,t);imshow(Dx(:,:,slice_show),[]); pause(0.2);
    title(['MRF-ICM����' num2str(t) '��' ]); 
    imwrite(Dx(:,:, slice_show) * 255, [path_folder_process, sprintf('image_4_ICM%d.png', t)]);
    imwrite(max(Dx,[],3) * 255, [path_folder_process, sprintf('image_4_ICM%d_projection.png', t)]);
    
    Dx_init = Dx;

    % ����ָ���̽��
    label_predict = zeros(size(label_save.img), 'uint8');
    label_predict(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2))...
        = uint8(Dx);
    label_save.img = label_predict;
    save_nii(label_save, [path_folder_process, 'GMM', MLE_org.distribution, 'ICM', num2str(t), '.nii.gz']);
    
    %% ����ICM���
    if do_crit
        fprintf('criterion ICM%d: ', t)
        eval(sprintf(...
            '[crit.ICM%d.accuracy, crit.ICM%d.sensitivity, crit.ICM%d.specificity, crit.ICM%d.precision, crit.ICM%d.dice, crit.ICM%d.loU] =  criterion(Dx, label_fg);',...
            t, t, t, t, t, t));
%         criterion(Dx, label_fg);
    end
    
end
% ����fig3
saveas(fig14, [path_folder_process, 'figure314, ICM-MRF process of post-segment bm4d'], 'fig');

%% ����ȥ��С��ͨ�򡣡�
Dx_removeisland = RemoveIsland(Dx, 10);
fprintf('criterion ICM remove island: ')
[crit.ICM_removeisland.accuracy, crit.ICM_removeisland.sensitivity, crit.ICM_removeisland.specificity, crit.ICM_removeisland.precision, ...
	crit.ICM_removeisland.dice, crit.ICM_removeisland.loU] = criterion(Dx_removeisland, label_fg);
% ����ͼƬ
imwrite(Dx_removeisland(:,:, slice_show) * 255, [path_folder_process, 'image_4_ICM_removeisland.png']);
imwrite(max(Dx_removeisland,[],3) * 255, [path_folder_process, 'image_4_ICM_removeisland_projection.png']);
% ����ָ���̽��
label_predict = zeros(size(label_save.img), 'uint8');
label_predict(image_box(1,1):image_box(1,2), image_box(2,1):image_box(2,2), image_box(3,1):image_box(3,2))...
    = uint8(Dx_removeisland);
label_save.img = label_predict;
save_nii(label_save, [path_folder_process, 'GMM', MLE_org.distribution, 'ICM_removeisland', num2str(t), '.nii.gz']);

%% ��ϸѪ�ܷ������
[label_fg_large,label_fg_small] = SeparationVesselSize(label_fg, 1);
[Dx_large, Dx_small] = SeparationVesselSize(Dx, 1);
[crit.ICM_large.accuracy, crit.ICM_large.sensitivity, crit.ICM_large.specificity, crit.ICM_large.precision, ...
	crit.ICM_large.dice, crit.ICM_large.loU] = criterion(Dx_large, label_fg_large);
[crit.ICM_small.accuracy, crit.ICM_small.sensitivity, crit.ICM_small.specificity, crit.ICM_small.precision, ...
	crit.ICM_small.dice, crit.ICM_small.loU] = criterion(Dx_small, label_fg_small);
[crit.ICM_large2orglabel.accuracy, crit.ICM_large2orglabel.sensitivity, crit.ICM_large2orglabel.specificity, crit.ICM_large2orglabel.precision, ...
	crit.ICM_large2orglabel.dice, crit.ICM_large2orglabel.loU] = criterion(Dx_large, label_fg);
% ����ͼƬ
imwrite(Dx_large(:,:, slice_show) * 255, [path_folder_process, 'image_5_ICM_large.png']);
imwrite(max(Dx_large,[],3) * 255, [path_folder_process, 'image_5_ICM_large_projection.png']);
imwrite(Dx_small(:,:, slice_show) * 255, [path_folder_process, 'image_5_ICM_small.png']);
imwrite(max(Dx_small,[],3) * 255, [path_folder_process, 'image_5_ICM_small_projection.png']);
imwrite(label_fg_large(:,:, slice_show) * 255, [path_folder_process, 'image_5_label_large.png']);
imwrite(max(label_fg_large,[],3) * 255, [path_folder_process, 'image_5_label_large_projection.png']);
imwrite(label_fg_small(:,:, slice_show) * 255, [path_folder_process, 'image_5_label_small.png']);
imwrite(max(label_fg_small,[],3) * 255, [path_folder_process, 'image_5_label_small_projection.png']);

%% ����record
record.parameters = parameters;
record.criterion = crit;
end