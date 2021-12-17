function [Dnew, image_pV, image_pB] = ICM_estimation(VBN,pxl_k,inPLX,Object,D,beta,Vx1,Vy1,Vz1,Iout,e,fg_ratio)
% VBN: EM�㷨�õ��Ĳ���
% pxl_k: �����ֲ��ĵ�ˮƽģ�͸���
% inPLX: ��ˮƽģ�͵õ������ر��
% Object: ǰ�����ڷֲ��еı��
% D: ����עͼ
% beta: 
% NB: 
% Vx1, Vy1, Vz1:
% Iout: ��ͨ��Ѫ���˲������
% e: �����������
W = VBN(:,3);
Li = 1:Object-1;%��Ŀ��������
sizeD = size(D);
Dnew = zeros(sizeD);
image_pV = zeros(sizeD);
image_pB = zeros(sizeD);
[a,b,c] = ind2sub(sizeD,inPLX);
s = [a b c];
 
% �����ˮƽģ�͸��ʣ�����������beta
pVy = pxl_k(:,Object)./W(Object);
pBy = sum(pxl_k(:,Li'),2)./sum(W(Li'));

if e == 2006 || e == 2006.1
    % [pB, pV] = clique_MPN(Object,s,D,beta,Vx1,Vy1,Vz1,Iout,e);
    [pB, pV] = clique_MPN(Object,inPLX,D,beta,Vx1,Vy1,Vz1,Iout,e,pVy,pBy,fg_ratio);
elseif e == 2013 || e == 2020 || e == 2013.1
    [pB, pV] = clique_MPN_zhou(Object,inPLX,D,beta,Vx1,Vy1,Vz1,Iout,e,pVy,pBy,fg_ratio);
end
% pB = pB(inPLX);
% pV = pV(inPLX);

post_V = pV.*pxl_k(:,Object)./W(Object);
post_B = pB.*sum(pxl_k(:,Li'),2)./sum(W(Li'));
Dnew(inPLX) = Object .* (post_V >post_B);
% Dnew(inPLX) = Object .* (pV >pB./2); % ����������EM�㷨�Ľ��

image_pV(inPLX) = pV;
image_pB(inPLX) = pB;
end
%****************** SubFunction cliqueMPN *************************
function [pB,pV] = clique_MPN(K,inPLX,D,beta,Vx,Vy,Vz,Iout,e,pVy,pBy,fg_ratio)
% 2006 M. Sabry Hassouna�ķ���

% ��D��26���򣬴���N��27���NΪinPLX�ĳ���
D_nb26 = GetNeighbors(D, inPLX);
D_nb26_foreground = D_nb26 == K;
D_nb26_background = D_nb26 ~= K;

% pVy��pBy��26����
image_pVy = zeros(size(D)); image_pBy = zeros(size(D));
image_pVy(inPLX) = pVy; image_pBy(inPLX) = pBy;
pVy_nb26 = GetNeighbors(image_pVy, inPLX);
pBy_nb26 = GetNeighbors(image_pBy, inPLX);
clear image_pVy image_pBy

if e==2006
    Uvw = 26 - sum(D_nb26_foreground, 2);% ��һ������������ƽ��
    Ubw = 26 - sum(D_nb26_background, 2);
elseif e==2006.1
    Uvw = 26 - sum(D_nb26_foreground.*pVy_nb26, 2); % ��һ����������������ϸ�ڸ���
    Ubw = 26 - sum(D_nb26_background.*pBy_nb26, 2);
end

Uv = beta(1).*Uvw;
Ub = beta(2).*Ubw;

%% ���ݵ�ˮƽģ�ͷ���beta
% �����߼�����������������x��
% X = V ���� U_v(x) - U_b(x) < 1/\beta (logP(Y=y|X=V) - logP(Y=y|X=B)) % Ҳ����һ��beta�����
%       ���� \beta_v��U_v(x) - \beta_b��U_b(x) <  logP(Y=y|X=V) - logP(Y=y|X=B)% ����beta�����
% ����������Ҫ�أ�U_v(x)��U_b(x)��logP(Y=y|X=V) - logP(Y=y|X=B)
pVy(pVy < 0.001) = 0.001; pBy(pBy < 0.001) = 0.001;
amplification_factor = 1000/max(beta(1) .* Uvw - beta(2) .* Ubw - (log(pVy)-log(pBy))); % ֱ��ͼ����������ķŴ�����
hist_Uvw = GetDecomposedHist(Uvw.*amplification_factor, D(inPLX), ones(size(Uvw)));
hist_Ubw = GetDecomposedHist(Ubw.*amplification_factor, D(inPLX), ones(size(Uvw)));
% hist_thrP = GetDecomposedHist((log(pVy)-log(pBy)).*amplification_factor, D(inPLX), ones(size(Uvw))); % ��error�������˳��������������ֵ����
% hist_ResultSegByZero = GetDecomposedHist(...
%     (Uvw - Ubw - 1 ./ beta(1) .* (log(pVy)-log(pBy))) .* amplification_factor, ...
%     D(inPLX), ones(size(Uvw)));
hist_ResultSegByZero = GetDecomposedHist(...
    (beta(1) .* Uvw - beta(2) .* Ubw - (log(pVy)-log(pBy))) .* amplification_factor, ...
    D(inPLX), ones(size(Uvw)));
% list_image = {'ResultSegByZero', 'Uvw', 'Ubw', 'thrP'};
% list_image = {'ResultSegByZero'};
list_image = {'ResultSegByZero', 'Uvw', 'Ubw'};
set(0,'DefaultFigureVisible', 'off'); % ����Ҫ��ʾ����Ҫ��ʾ����Ҫ��ʾ��
figure(31); 
fig31 = gcf;
set(fig31,'Units','centimeter','Position',[5 5 25 10]);
for j = 1:length(list_image)
    eval(['hist_temp = hist_', char(list_image(j)), ';']);

    subplot(1,length(list_image),j); hold on
    title(strrep(char(list_image(j)), '_', ' '));
    % ��ǰ��/����/ȫ��ֱ��ͼ
    plot(hist_temp.X_mask / amplification_factor, hist_temp.N_mask / numel(inPLX)); hold on
    plot(hist_temp.X_fg / amplification_factor, hist_temp.N_fg / numel(inPLX));
    plot(hist_temp.X_bg / amplification_factor, hist_temp.N_bg / numel(inPLX));
    xlim([-10 10]);
    ylim([0 0.001]);
    legend(sprintf('liver: fg/bg = %.3f/%.3f', ...% ���˳�������ԭ���ı�ǩ�����ϣ��±�ǩ������ǰ��/�����ĸ���
        sum(hist_temp.N_mask(hist_temp.X_mask / amplification_factor < 0)) / numel(inPLX),...
        sum(hist_temp.N_mask(hist_temp.X_mask / amplification_factor > 0)) / numel(inPLX)), ...
            sprintf('vessel: fg/bg = %.3f/%.3f', ...
        sum(hist_temp.N_fg(hist_temp.X_fg / amplification_factor < 0)) / numel(inPLX),...
        sum(hist_temp.N_fg(hist_temp.X_fg / amplification_factor > 0)) / numel(inPLX)), ...
            sprintf('background: fg/bg = %.3f/%.3f', ...
        sum(hist_temp.N_bg(hist_temp.X_bg / amplification_factor < 0)) / numel(inPLX),...
        sum(hist_temp.N_bg(hist_temp.X_bg / amplification_factor > 0)) / numel(inPLX)));
    hold off
end

% ����У������C����
% X = V ���� logP(Y=y|X=V) - logP(Y=y|X=B) - \beta_v��U_v(x) + \beta_b��U_b(x) > 0
% �Ľ�Ϊ
% X = V ���� logP(Y=y|X=V) - logP(Y=y|X=B) - \beta_v��U_v(x) + \beta_b��U_b(x) > C
% ʹP(X=V)=ǰ���ı���fg_ratio
[~,idx] = min(abs(cumsum(hist_ResultSegByZero.N_mask) / numel(inPLX) - fg_ratio));
C = hist_ResultSegByZero.X_mask(idx) / amplification_factor;

% pV = exp(-Uv);
% pB = exp(-Ub);
pV = exp(-Uv+C);
pB = exp(-Ub);
end
%****************** SubFunction cliqueMPN *************************
% function [pB,pV,NsigmaV,NsigmaB] = clique_MPN(K,s,D,beta,Iout,Vx,Vy,Vz,e)
function [pB,pV,NsigmaV,NsigmaB] = clique_MPN_zhou(K,inPLX,D,beta,Vx,Vy,Vz,Iout,e,pVy,pBy,fg_ratio)

% �ڼȶ���26��������D_nb26�У�����ռ�����ľ���ṹ
A = [5 11 13 15 17 23]; % D_nb26�е�6�����㣻
% ��(i,j,n)Ϊ���ģ���26������������ɺ���ǰ���������ҡ��������Ϸ�Ϊ8�������壬ÿ���������(i,j,n)����7������,������Ķ�����ֵö�����£�
C = [1 2 4 5 10 11 13;2 3 5 6 11 12 15;4 5 7 8 13 16 17;5 6 8 9 15 17 18;...
     10 11 13 19 20 22 23;11 12 15 20 21 23 24;13 16 17 22 23 25 26;15 17 18 23 24 26 27];
% ��(i,j,n)Ϊ���ģ���26������������ɺ���ǰ���������ҡ��������Ϸ�Ϊ6�������壬ÿ���������(i,j,n)����9������,������Ķ�����ֵö�����£�
F = [4 5 10 11 13 16 17 22 23;5 6 11 12 15 17 18 23 24;2 5 10 11 12 13 15 20 23; ...
     5 8 13 15 16 17 18 23 26;2 4 5 6 8 11 13 15 17;11 13 15 17 20 22 23 24 26]; 

% ��D��26���򣬴���N��27���NΪinPLX�ĳ���
D_nb26 = GetNeighbors(D, inPLX);
D_nb26_foreground = D_nb26 == K;
D_nb26_background = D_nb26 ~= K;

if e == 2020
    % Vx, Vy, Vz, Iout�������N��27
    Vx_nb26 = GetNeighbors(Vx, inPLX);
    Vy_nb26 = GetNeighbors(Vy, inPLX);
    Vz_nb26 = GetNeighbors(Vz, inPLX);
    Iout_nb26 = GetNeighbors(Iout, inPLX);
    % Vx, Vy, Vz������ɵ���������СΪN��3��27
    V = zeros(numel(inPLX), 3, 27); 
    V(:,1,:) = Vx_nb26; V(:,2,:) = Vy_nb26; V(:,3,:) = Vz_nb26;
    % Vx, Vy, Vz��ɵ���������СΪN��3��1��Ϊ�˷����ǰ���ˣ���չΪN��3��27
    Vm = V(:,:,14);
    Vm = repmat(Vm, 1, 1, 27);
elseif e == 2013.1
    % pVy��pBy�������N��27
    image_pVy = zeros(size(D)); image_pBy = zeros(size(D));
    image_pVy(inPLX) = pVy; image_pBy(inPLX) = pBy;
    pVy_nb26 = GetNeighbors(image_pVy, inPLX);
    pBy_nb26 = GetNeighbors(image_pBy, inPLX);
    clear image_pVy image_pBy
end

if e == 2013 % 2013�汾
    % 26������ά�ռ��N�����ż��ϣ����ﲻֻ�������ˣ�����inPLX���ȵı�
    %A
    NsigmaV_6 = sum(D_nb26_foreground(:,A), 2);
    NsigmaB_6 = sum(D_nb26_background(:,A), 2);
    %C
    NsigmaV_7 = zeros(length(inPLX), 1);
    NsigmaB_7 = zeros(length(inPLX), 1);
    for i = 1:8
        NsigmaV_7 = max(sum(D_nb26_foreground(:,squeeze(C(i,:))), 2), NsigmaV_7);
        NsigmaB_7 = max(sum(D_nb26_background(:,squeeze(C(i,:))), 2), NsigmaB_7);
    end
    %F
    NsigmaV_9 = zeros(length(inPLX), 1);
    NsigmaB_9 = zeros(length(inPLX), 1);
    for i = 1:6
        NsigmaV_9 = max(sum(D_nb26_foreground(:,squeeze(F(i,:))), 2), NsigmaV_9);
        NsigmaB_9 = max(sum(D_nb26_background(:,squeeze(F(i,:))), 2), NsigmaB_9);
    end
elseif e == 2020 % 2020�汾
%     epsilon1 = 0.5; epsilon2 = 0.5;
%     epsilon1 = 1; epsilon2 = 0;
    epsilon1 = 0; epsilon2 = 1;
    costheta = sum(Vm .* V, 2); % ���ĵĵ����������ļнǵ�cos����СΪN��27
    dI = abs(Iout_nb26 - Iout_nb26(:,14)); % ���ĵĵ�����������|��I| ���˲�ǿ�Ȳ��СΪN��27
    % A
%     costheta_6 = costheta(:,A); % ȡ��״A�Ĳ��ּ���
%     S = 1- abs(costheta_6); % similarity�������Ƕȼнǡ���costheta����S��
%     e1 = exp(S); % �ҷ���������û�ⲽ�����ҡ��ⲻ��ѧ������S�ķ�Χ����1-e��������0-1��ֱ�ӵ���rho��ֵ�϶�����1����ô1-rho��С��0�ˡ�
%     dV = abs(Iout_nb26 - Iout_nb26(:,14)); % |��I| ��ǿ�Ȳ����dV��
%     dV_norm = dV ./ (max(dV, [], 2) + eps); % |��I| / (|��Imax|)
%     rho = epsilon1 * dV_norm + epsilon2 * e1; % �нǡ���ǿ�Ȳ����rho��
%     NsigmaV_6 = sum(D_nb26_background(:,A) .* rho, 2); % xs, xs'ͬΪL_V��0��xs, xs'��ͬ���rho
%     NsigmaB_6 = sum(D_nb26_foreground(:,A) .* rho + D_nb26_background(:,A) .* abs(1 - rho), 2); %��1 - rho�ķ�Χ��������0-1������ô����������
    [NsigmaV_6, NsigmaB_6] = GetNsigma2020(costheta, dI, A, D_nb26_foreground, D_nb26_background, epsilon1, epsilon2);
    % C
    [NsigmaV_7, NsigmaB_7] = GetNsigma2020(costheta, dI, C, D_nb26_foreground, D_nb26_background, epsilon1, epsilon2);
    % F
    [NsigmaV_9, NsigmaB_9] = GetNsigma2020(costheta, dI, F, D_nb26_foreground, D_nb26_background, epsilon1, epsilon2);
elseif e == 2013.1
    % A
    [NsigmaV_6, NsigmaB_6] = GetNsigma2021(pVy_nb26, pBy_nb26, A, D_nb26_foreground, D_nb26_background);
    % C
    [NsigmaV_7, NsigmaB_7] = GetNsigma2021(pVy_nb26, pBy_nb26, C, D_nb26_foreground, D_nb26_background);
    % F
    [NsigmaV_9, NsigmaB_9] = GetNsigma2021(pVy_nb26, pBy_nb26, F, D_nb26_foreground, D_nb26_background);

end
if e == 2013
    % D_nb26�о���A���������Ŀ��ͱ�������
    Uv6 = 6 - NsigmaV_6; Ub6 = 6 - NsigmaB_6;
    % D_nb26�о���C���������Ŀ��ͱ�������
    Uv_bound1 = 7 - NsigmaV_7; Ub_bound1 = 7 - NsigmaB_7;
    Uv_bound2 = 9 - NsigmaV_9; Ub_bound2 = 9 - NsigmaB_9;
elseif e == 2020 || 2013.1
    % D_nb26�о���A���������Ŀ��ͱ�������
    Uv6 = NsigmaV_6; Ub6 = NsigmaB_6;
    % D_nb26�о���C���������Ŀ��ͱ�������
    Uv_bound1 = NsigmaV_7; Ub_bound1 = NsigmaB_7;
    Uv_bound2 = NsigmaV_9; Ub_bound2 = NsigmaB_9;
end
% ����Ŀ����ʣ�����D(i,j,n)�Ƿ�ΪK����Χ���K�϶�ʱ��������С��Ŀ��������
% Uvw = min([Uv6 Uv_bound1 Uv_bound2]);%Uvw = [Uv6 Uv_bound1 Uv_bound2];
Uvw = min(cat(2, Uv6, Uv_bound1, Uv_bound2), [], 2);%Uvw = [Uv6 Uv_bound1 Uv_bound2];
Uv = beta(1).*Uvw;
% ���㱳�����ʣ�����Χ��ǲ�ΪK�Ľ϶�ʱ��������С�������������
Ubw = min(cat(2, Ub6, Ub_bound1, Ub_bound2), [], 2);%Ubw = [Ub6 Ub_bound1 Ub_bound2];
Ub = beta(2).*Ubw;

%% ���ݵ�ˮƽģ�ͷ���beta
% �����߼�����������������x��
% X = V ���� U_v(x) - U_b(x) < 1/\beta (logP(Y=y|X=V) - logP(Y=y|X=B)) % Ҳ����һ��beta�����
%       ���� \beta_v��U_v(x) - \beta_b��U_b(x) <  logP(Y=y|X=V) - logP(Y=y|X=B)% ����beta�����
% ����������Ҫ�أ�U_v(x)��U_b(x)��logP(Y=y|X=V) - logP(Y=y|X=B)
pVy(pVy < 0.001) = 0.001; pBy(pBy < 0.001) = 0.001;
amplification_factor = 1000/max(beta(1) .* Uvw - beta(2) .* Ubw - (log(pVy)-log(pBy))); % ֱ��ͼ����������ķŴ�����
hist_Uvw = GetDecomposedHist(Uvw.*amplification_factor, D(inPLX), ones(size(Uvw)));
hist_Ubw = GetDecomposedHist(Ubw.*amplification_factor, D(inPLX), ones(size(Uvw)));
% hist_thrP = GetDecomposedHist((log(pVy)-log(pBy)).*amplification_factor, D(inPLX), ones(size(Uvw))); % ��error�������˳��������������ֵ����
% hist_ResultSegByZero = GetDecomposedHist(...
%     (Uvw - Ubw - 1 ./ beta(1) .* (log(pVy)-log(pBy))) .* amplification_factor, ...
%     D(inPLX), ones(size(Uvw)));
hist_ResultSegByZero = GetDecomposedHist(...
    (beta(1) .* Uvw - beta(2) .* Ubw - (log(pVy)-log(pBy))) .* amplification_factor, ...
    D(inPLX), ones(size(Uvw)));
% list_image = {'ResultSegByZero', 'Uvw', 'Ubw', 'thrP'};
% list_image = {'ResultSegByZero'};
list_image = {'ResultSegByZero', 'Uvw', 'Ubw'};
figure(31)
fig31 = gcf;
set(fig31,'Units','centimeter','Position',[5 5 25 10]);
for j = 1:length(list_image)
    eval(['hist_temp = hist_', char(list_image(j)), ';']);

    subplot(1,length(list_image),j); hold on
    title(strrep(char(list_image(j)), '_', ' '));
    % ��ǰ��/����/ȫ��ֱ��ͼ
    plot(hist_temp.X_mask / amplification_factor, hist_temp.N_mask / numel(inPLX)); hold on
    plot(hist_temp.X_fg / amplification_factor, hist_temp.N_fg / numel(inPLX));
    plot(hist_temp.X_bg / amplification_factor, hist_temp.N_bg / numel(inPLX));
    xlim([-10 10]);
    ylim([0 0.001]);
    legend(sprintf('liver: fg/bg = %.3f/%.3f', ...% ���˳�������ԭ���ı�ǩ�����ϣ��±�ǩ������ǰ��/�����ĸ���
        sum(hist_temp.N_mask(hist_temp.X_mask / amplification_factor < 0)) / numel(inPLX),...
        sum(hist_temp.N_mask(hist_temp.X_mask / amplification_factor > 0)) / numel(inPLX)), ...
            sprintf('vessel: fg/bg = %.3f/%.3f', ...
        sum(hist_temp.N_fg(hist_temp.X_fg / amplification_factor < 0)) / numel(inPLX),...
        sum(hist_temp.N_fg(hist_temp.X_fg / amplification_factor > 0)) / numel(inPLX)), ...
            sprintf('background: fg/bg = %.3f/%.3f', ...
        sum(hist_temp.N_bg(hist_temp.X_bg / amplification_factor < 0)) / numel(inPLX),...
        sum(hist_temp.N_bg(hist_temp.X_bg / amplification_factor > 0)) / numel(inPLX)));
    hold off
end

% ����У������C����
% X = V ���� logP(Y=y|X=V) - logP(Y=y|X=B) - \beta_v��U_v(x) + \beta_b��U_b(x) > 0
% �Ľ�Ϊ
% X = V ���� logP(Y=y|X=V) - logP(Y=y|X=B) - \beta_v��U_v(x) + \beta_b��U_b(x) > C
% ʹP(X=V)=ǰ���ı���fg_ratio
[~,idx] = min(abs(cumsum(hist_ResultSegByZero.N_mask) / numel(inPLX) - fg_ratio));
C = hist_ResultSegByZero.X_mask(idx) / amplification_factor;

% pV = exp(-Uv);
% pB = exp(-Ub);
pV = exp(-Uv+C);
pB = exp(-Ub);
% ����Ŀ�����NsigmaV�ͱ�������NsigmaB
% NsigmaV = min([NsigmaV_6 NsigmaV_7 NsigmaV_9]);
% NsigmaB = min([NsigmaB_6 NsigmaB_7 NsigmaB_9]);
NsigmaV = min(cat(2, NsigmaV_6, NsigmaV_7, NsigmaV_9), [], 2);
NsigmaB = min(cat(2, NsigmaB_6, NsigmaB_7, NsigmaB_9), [], 2);
end

function D_nb26 = GetNeighbors(D, inPLX)
% ��D��26���򣬴���N��27���NΪD����������
sizeD = size(D);
D_expand = zeros(size(D)+2);
D_expand(2:sizeD(1)+1, 2:sizeD(2)+1, 2:sizeD(3)+1) = D;
D_neighbors = zeros([sizeD(1), sizeD(2), sizeD(3), 27]);
for i = [-1, 0, 1]
    for j = [-1, 0, 1]
        for k = [-1, 0, 1]
            D_neighbors(:,:,:, (i+1) + (j+1)*3 + (k+1)*9 + 1) = D_expand((2:sizeD(1)+1)+i, (2:sizeD(2)+1)+j, (2:sizeD(3)+1)+k);
        end
    end
end
D_neighbors = reshape(D_neighbors, numel(D), 27);
% ֻ����inPLXָʾ�Ĳ���
D_nb26 = D_neighbors(inPLX, :);
end

function [NsigmaV, NsigmaB] = GetNsigma2020(costheta, dI, NBS, D_nb26_foreground, D_nb26_background, epsilon1, epsilon2)
% 2020 lina�汾��Nsigma�󷨣���\scriptP_C_(m,j)(x_s, x_s')
% INPUT
% costheta��N��27����ÿ����������26�������صļнǵ�cos��PS: �����������м䣬��14����
% Iout_nb26��N��27����ÿ����������26�������صĶ�ͨ���˲����
% NBS��NBS��ת��������NBS�����������󣬰������NBS���������ת����
% D_nb26_foreground��N��27����ÿ����������26��������Ϊǰ�����Ϊ1������Ϊ0
% D_nb26_background��N��27����ÿ����������26��������Ϊ�������Ϊ1������Ϊ0
% epsilon1, epsilon2����Ϊ1����������PWS-potential����Ԫ��|��I| / (|��Imax|)��Similarity��Ȩ

[NBS_num, ~] = size(NBS); % ���class NBS������
[lenPLX, ~] = size(costheta); % ������������
% NsigmaV = zeros(lenPLX, 1);
% NsigmaB = zeros(lenPLX, 1);
nsigmaV = zeros(lenPLX, NBS_num);
nsigmaB = zeros(lenPLX, NBS_num);
for j = 1:NBS_num
    costheta_NBSj = costheta(:,NBS(j,:));  % ȡ��״NBS�Ĳ��ּ���
    S = 1- abs(costheta_NBSj); % similarity�������Ƕȼнǡ���costheta����S��
    e1 = exp(S); % �ҷ���������û�ⲽ�����ҡ��ⲻ��ѧ������S�ķ�Χ����1-e��������0-1��ֱ�ӵ���rho��ֵ�϶�����1����ô1-rho��С��0�ˡ�
    
    dI_NBSj = dI(:,NBS(j,:)); % |��I| ��ǿ�Ȳ����dV��
    dV_norm = dI_NBSj ./ (max(dI_NBSj, [], 2) + eps);  % |��I| / (|��Imax|)
    % ԭ����
    rho = epsilon1 * dV_norm + epsilon2 * e1; % �нǡ���ǿ�Ȳ����rho��
    
%     NsigmaV = min(nsigmaV, sum(D_nb26_background(:,NBS(j,:)) .* rho, 2)); % xs, xs'ͬΪL_V��0��xs, xs'��ͬ���rho
%     NsigmaB = min(nsigmaB, sum(D_nb26_foreground(:,NBS(j,:)) .* rho + D_nb26_background(:,NBS(j,:)) .* abs(1 - rho), 2)); %��1 - rho�ķ�Χ��������0-1������ô����������
    nsigmaV(:,j) = sum(D_nb26_background(:,NBS(j,:)) .* rho, 2); % xs, xs'ͬΪL_V��0��xs, xs'��ͬ���rho
    nsigmaB(:,j) = sum(D_nb26_foreground(:,NBS(j,:)) .* rho + D_nb26_background(:,NBS(j,:)) .* abs(1 - rho), 2); %��1 - rho�ķ�Χ��������0-1������ô����������
    % ԭ����
%     % �¡����ｻ����rho��similarity��|��I|��Ȩ�أ�
%     rho = 0 * dV_norm + 1 * e1; % �нǡ���ǿ�Ȳ����rho��
%     nsigmaV(:,j) = sum(D_nb26_background(:,NBS(j,:)) .* rho, 2); % xs, xs'ͬΪL_V��0��xs, xs'��ͬ���rho
%     rho = 0.5 * dV_norm + 0.5 * e1; % �нǡ���ǿ�Ȳ����rho��
%     nsigmaB(:,j) = sum(D_nb26_foreground(:,NBS(j,:)) .* rho + D_nb26_background(:,NBS(j,:)) .* abs(1 - rho), 2); %��1 - rho�ķ�Χ��������0-1������ô����������
%     % �¡�
end
NsigmaV = min(nsigmaV, [], 2);
NsigmaB = min(nsigmaB, [], 2);
end

function [NsigmaV, NsigmaB] = GetNsigma2021(pVy, pBy, NBS, D_nb26_foreground, D_nb26_background)
% 2020 lina�汾��Nsigma�󷨣���\scriptP_C_(m,j)(x_s, x_s')
% INPUT
% pVy��N��27����ÿ����������26�������ص�pVy��PS: �����������м䣬��14����
% pBy��N��27����ÿ����������26�������ص�pBy
% NBS��NBS��ת��������NBS�����������󣬰������NBS���������ת����
% D_nb26_foreground��N��27����ÿ����������26��������Ϊǰ�����Ϊ1������Ϊ0
% D_nb26_background��N��27����ÿ����������26��������Ϊ�������Ϊ1������Ϊ0
% epsilon1, epsilon2����Ϊ1����������PWS-potential����Ԫ��|��I| / (|��Imax|)��Similarity��Ȩ

[NBS_num, ~] = size(NBS); % ���class NBS������
[lenPLX, ~] = size(pVy); % ������������
% NsigmaV = zeros(lenPLX, 1);
% NsigmaB = zeros(lenPLX, 1);
nsigmaV = zeros(lenPLX, NBS_num);
nsigmaB = zeros(lenPLX, NBS_num);
for j = 1:NBS_num
    pV = pVy(:,NBS(j,:)) ./ (pVy(:,NBS(j,:)) + pBy(:,NBS(j,:)));
    % ԭ����
    rho = pV; % pV����rho��
   nsigmaV(:,j) = sum(D_nb26_background(:,NBS(j,:)) .* rho, 2); % xs, xs'ͬΪL_V��0��xs, xs'��ͬ���rho
    nsigmaB(:,j) = sum(D_nb26_foreground(:,NBS(j,:)) .* rho + D_nb26_background(:,NBS(j,:)) .* abs(1 - rho), 2); %��1 - rho�ķ�Χ��������0-1������ô����������
end
NsigmaV = min(nsigmaV, [], 2);
NsigmaB = min(nsigmaB, [], 2);
end