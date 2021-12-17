function [Dnew, image_pV, image_pB] = ICM_estimation(VBN,pxl_k,inPLX,Object,D,beta,Vx1,Vy1,Vz1,Iout,e,fg_ratio)
% VBN: EM算法得到的参数
% pxl_k: 各个分布的低水平模型概率
% inPLX: 低水平模型得到的体素编号
% Object: 前景类在分布中的编号
% D: 类别标注图
% beta: 
% NB: 
% Vx1, Vy1, Vz1:
% Iout: 多通道血管滤波器结果
% e: 能量函数编号
W = VBN(:,3);
Li = 1:Object-1;%非目标类的序号
sizeD = size(D);
Dnew = zeros(sizeD);
image_pV = zeros(sizeD);
image_pB = zeros(sizeD);
[a,b,c] = ind2sub(sizeD,inPLX);
s = [a b c];
 
% 计算低水平模型概率，以用来估计beta
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
% Dnew(inPLX) = Object .* (pV >pB./2); % 这里无视了EM算法的结果

image_pV(inPLX) = pV;
image_pB(inPLX) = pB;
end
%****************** SubFunction cliqueMPN *************************
function [pB,pV] = clique_MPN(K,inPLX,D,beta,Vx,Vy,Vz,Iout,e,pVy,pBy,fg_ratio)
% 2006 M. Sabry Hassouna的方法

% 找D的26邻域，存在N×27表里，N为inPLX的长度
D_nb26 = GetNeighbors(D, inPLX);
D_nb26_foreground = D_nb26 == K;
D_nb26_background = D_nb26 ~= K;

% pVy和pBy的26邻域
image_pVy = zeros(size(D)); image_pBy = zeros(size(D));
image_pVy(inPLX) = pVy; image_pBy(inPLX) = pBy;
pVy_nb26 = GetNeighbors(image_pVy, inPLX);
pBy_nb26 = GetNeighbors(image_pBy, inPLX);
clear image_pVy image_pBy

if e==2006
    Uvw = 26 - sum(D_nb26_foreground, 2);% 这一种能量函数更平滑
    Ubw = 26 - sum(D_nb26_background, 2);
elseif e==2006.1
    Uvw = 26 - sum(D_nb26_foreground.*pVy_nb26, 2); % 这一种能量函数保留的细节更多
    Ubw = 26 - sum(D_nb26_background.*pBy_nb26, 2);
end

Uv = beta(1).*Uvw;
Ub = beta(2).*Ubw;

%% 根据低水平模型分析beta
% 根据逻辑推理，对于所有像素x，
% X = V ←→ U_v(x) - U_b(x) < 1/\beta (logP(Y=y|X=V) - logP(Y=y|X=B)) % 也就是一个beta的情况
%       ←→ \beta_v·U_v(x) - \beta_b·U_b(x) <  logP(Y=y|X=V) - logP(Y=y|X=B)% 两个beta的情况
% 分析这三个要素：U_v(x)，U_b(x)，logP(Y=y|X=V) - logP(Y=y|X=B)
pVy(pVy < 0.001) = 0.001; pBy(pBy < 0.001) = 0.001;
amplification_factor = 1000/max(beta(1) .* Uvw - beta(2) .* Ubw - (log(pVy)-log(pBy))); % 直方图分析横坐标的放大因子
hist_Uvw = GetDecomposedHist(Uvw.*amplification_factor, D(inPLX), ones(size(Uvw)));
hist_Ubw = GetDecomposedHist(Ubw.*amplification_factor, D(inPLX), ones(size(Uvw)));
% hist_thrP = GetDecomposedHist((log(pVy)-log(pBy)).*amplification_factor, D(inPLX), ones(size(Uvw))); % 【error：超出了程序允许的最大变量值。】
% hist_ResultSegByZero = GetDecomposedHist(...
%     (Uvw - Ubw - 1 ./ beta(1) .* (log(pVy)-log(pBy))) .* amplification_factor, ...
%     D(inPLX), ones(size(Uvw)));
hist_ResultSegByZero = GetDecomposedHist(...
    (beta(1) .* Uvw - beta(2) .* Ubw - (log(pVy)-log(pBy))) .* amplification_factor, ...
    D(inPLX), ones(size(Uvw)));
% list_image = {'ResultSegByZero', 'Uvw', 'Ubw', 'thrP'};
% list_image = {'ResultSegByZero'};
list_image = {'ResultSegByZero', 'Uvw', 'Ubw'};
set(0,'DefaultFigureVisible', 'off'); % 【不要显示，不要显示，不要显示】
figure(31); 
fig31 = gcf;
set(fig31,'Units','centimeter','Position',[5 5 25 10]);
for j = 1:length(list_image)
    eval(['hist_temp = hist_', char(list_image(j)), ';']);

    subplot(1,length(list_image),j); hold on
    title(strrep(char(list_image(j)), '_', ' '));
    % 画前景/背景/全部直方图
    plot(hist_temp.X_mask / amplification_factor, hist_temp.N_mask / numel(inPLX)); hold on
    plot(hist_temp.X_fg / amplification_factor, hist_temp.N_fg / numel(inPLX));
    plot(hist_temp.X_bg / amplification_factor, hist_temp.N_bg / numel(inPLX));
    xlim([-10 10]);
    ylim([0 0.001]);
    legend(sprintf('liver: fg/bg = %.3f/%.3f', ...% 这儿顺便计算在原来的标签基础上，新标签里面是前景/背景的概率
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

% 计算校正参数C，将
% X = V ←→ logP(Y=y|X=V) - logP(Y=y|X=B) - \beta_v·U_v(x) + \beta_b·U_b(x) > 0
% 改进为
% X = V ←→ logP(Y=y|X=V) - logP(Y=y|X=B) - \beta_v·U_v(x) + \beta_b·U_b(x) > C
% 使P(X=V)=前景的比例fg_ratio
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

% 在既定的26邻域数组D_nb26中，定义空间邻域的矩阵结构
A = [5 11 13 15 17 23]; % D_nb26中的6邻域格点；
% 以(i,j,n)为中心，将26邻域的立方体由后向前、由左至右、由下至上分为8个立方体，每个立方体除(i,j,n)外有7个顶点,立方体的顶点标记值枚举如下：
C = [1 2 4 5 10 11 13;2 3 5 6 11 12 15;4 5 7 8 13 16 17;5 6 8 9 15 17 18;...
     10 11 13 19 20 22 23;11 12 15 20 21 23 24;13 16 17 22 23 25 26;15 17 18 23 24 26 27];
% 以(i,j,n)为中心，将26邻域的立方体由后向前、由左至右、由下至上分为6个六面体，每个六面体除(i,j,n)外有9个顶点,六面体的顶点标记值枚举如下：
F = [4 5 10 11 13 16 17 22 23;5 6 11 12 15 17 18 23 24;2 5 10 11 12 13 15 20 23; ...
     5 8 13 15 16 17 18 23 26;2 4 5 6 8 11 13 15 17;11 13 15 17 20 22 23 24 26]; 

% 找D的26邻域，存在N×27表里，N为inPLX的长度
D_nb26 = GetNeighbors(D, inPLX);
D_nb26_foreground = D_nb26 == K;
D_nb26_background = D_nb26 ~= K;

if e == 2020
    % Vx, Vy, Vz, Iout的邻域表N×27
    Vx_nb26 = GetNeighbors(Vx, inPLX);
    Vy_nb26 = GetNeighbors(Vy, inPLX);
    Vz_nb26 = GetNeighbors(Vz, inPLX);
    Iout_nb26 = GetNeighbors(Iout, inPLX);
    % Vx, Vy, Vz邻域组成的向量表，大小为N×3×27
    V = zeros(numel(inPLX), 3, 27); 
    V(:,1,:) = Vx_nb26; V(:,2,:) = Vy_nb26; V(:,3,:) = Vz_nb26;
    % Vx, Vy, Vz组成的向量表，大小为N×3×1，为了方便和前面点乘，扩展为N×3×27
    Vm = V(:,:,14);
    Vm = repmat(Vm, 1, 1, 27);
elseif e == 2013.1
    % pVy和pBy的邻域表N×27
    image_pVy = zeros(size(D)); image_pBy = zeros(size(D));
    image_pVy(inPLX) = pVy; image_pBy(inPLX) = pBy;
    pVy_nb26 = GetNeighbors(image_pVy, inPLX);
    pBy_nb26 = GetNeighbors(image_pBy, inPLX);
    clear image_pVy image_pBy
end

if e == 2013 % 2013版本
    % 26邻域三维空间的N阶势团集合，这里不只是数字了，而是inPLX长度的表
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
elseif e == 2020 % 2020版本
%     epsilon1 = 0.5; epsilon2 = 0.5;
%     epsilon1 = 1; epsilon2 = 0;
    epsilon1 = 0; epsilon2 = 1;
    costheta = sum(Vm .* V, 2); % 中心的点与它邻域点的夹角的cos，大小为N×27
    dI = abs(Iout_nb26 - Iout_nb26(:,14)); % 中心的点与它邻域点的|▽I| ，滤波强度差，大小为N×27
    % A
%     costheta_6 = costheta(:,A); % 取形状A的部分计算
%     S = 1- abs(costheta_6); % similarity，两个角度夹角↓→costheta↑→S↓
%     e1 = exp(S); % 我发誓论文里没这步，而且【这不科学，这样S的范围就是1-e，而不是0-1，直接导致rho的值肯定大于1，那么1-rho就小于0了】
%     dV = abs(Iout_nb26 - Iout_nb26(:,14)); % |▽I| ，强度差↓，dV↓
%     dV_norm = dV ./ (max(dV, [], 2) + eps); % |▽I| / (|▽Imax|)
%     rho = epsilon1 * dV_norm + epsilon2 * e1; % 夹角↓，强度差↓，rho↓
%     NsigmaV_6 = sum(D_nb26_background(:,A) .* rho, 2); % xs, xs'同为L_V→0，xs, xs'不同类→rho
%     NsigmaB_6 = sum(D_nb26_foreground(:,A) .* rho + D_nb26_background(:,A) .* abs(1 - rho), 2); %【1 - rho的范围根本不是0-1，这怎么能这样！】
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
    % D_nb26中矩阵A表达的邻域的目标和背景能量
    Uv6 = 6 - NsigmaV_6; Ub6 = 6 - NsigmaB_6;
    % D_nb26中矩阵C表达的邻域的目标和背景能量
    Uv_bound1 = 7 - NsigmaV_7; Ub_bound1 = 7 - NsigmaB_7;
    Uv_bound2 = 9 - NsigmaV_9; Ub_bound2 = 9 - NsigmaB_9;
elseif e == 2020 || 2013.1
    % D_nb26中矩阵A表达的邻域的目标和背景能量
    Uv6 = NsigmaV_6; Ub6 = NsigmaB_6;
    % D_nb26中矩阵C表达的邻域的目标和背景能量
    Uv_bound1 = NsigmaV_7; Ub_bound1 = NsigmaB_7;
    Uv_bound2 = NsigmaV_9; Ub_bound2 = NsigmaB_9;
end
% 计算目标概率：无论D(i,j,n)是否为K，周围标记K较多时，能量最小，目标概率最大
% Uvw = min([Uv6 Uv_bound1 Uv_bound2]);%Uvw = [Uv6 Uv_bound1 Uv_bound2];
Uvw = min(cat(2, Uv6, Uv_bound1, Uv_bound2), [], 2);%Uvw = [Uv6 Uv_bound1 Uv_bound2];
Uv = beta(1).*Uvw;
% 计算背景概率：当周围标记不为K的较多时，能量最小，背景概率最大
Ubw = min(cat(2, Ub6, Ub_bound1, Ub_bound2), [], 2);%Ubw = [Ub6 Ub_bound1 Ub_bound2];
Ub = beta(2).*Ubw;

%% 根据低水平模型分析beta
% 根据逻辑推理，对于所有像素x，
% X = V ←→ U_v(x) - U_b(x) < 1/\beta (logP(Y=y|X=V) - logP(Y=y|X=B)) % 也就是一个beta的情况
%       ←→ \beta_v·U_v(x) - \beta_b·U_b(x) <  logP(Y=y|X=V) - logP(Y=y|X=B)% 两个beta的情况
% 分析这三个要素：U_v(x)，U_b(x)，logP(Y=y|X=V) - logP(Y=y|X=B)
pVy(pVy < 0.001) = 0.001; pBy(pBy < 0.001) = 0.001;
amplification_factor = 1000/max(beta(1) .* Uvw - beta(2) .* Ubw - (log(pVy)-log(pBy))); % 直方图分析横坐标的放大因子
hist_Uvw = GetDecomposedHist(Uvw.*amplification_factor, D(inPLX), ones(size(Uvw)));
hist_Ubw = GetDecomposedHist(Ubw.*amplification_factor, D(inPLX), ones(size(Uvw)));
% hist_thrP = GetDecomposedHist((log(pVy)-log(pBy)).*amplification_factor, D(inPLX), ones(size(Uvw))); % 【error：超出了程序允许的最大变量值。】
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
    % 画前景/背景/全部直方图
    plot(hist_temp.X_mask / amplification_factor, hist_temp.N_mask / numel(inPLX)); hold on
    plot(hist_temp.X_fg / amplification_factor, hist_temp.N_fg / numel(inPLX));
    plot(hist_temp.X_bg / amplification_factor, hist_temp.N_bg / numel(inPLX));
    xlim([-10 10]);
    ylim([0 0.001]);
    legend(sprintf('liver: fg/bg = %.3f/%.3f', ...% 这儿顺便计算在原来的标签基础上，新标签里面是前景/背景的概率
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

% 计算校正参数C，将
% X = V ←→ logP(Y=y|X=V) - logP(Y=y|X=B) - \beta_v·U_v(x) + \beta_b·U_b(x) > 0
% 改进为
% X = V ←→ logP(Y=y|X=V) - logP(Y=y|X=B) - \beta_v·U_v(x) + \beta_b·U_b(x) > C
% 使P(X=V)=前景的比例fg_ratio
[~,idx] = min(abs(cumsum(hist_ResultSegByZero.N_mask) / numel(inPLX) - fg_ratio));
C = hist_ResultSegByZero.X_mask(idx) / amplification_factor;

% pV = exp(-Uv);
% pB = exp(-Ub);
pV = exp(-Uv+C);
pB = exp(-Ub);
% 计算目标点数NsigmaV和背景点数NsigmaB
% NsigmaV = min([NsigmaV_6 NsigmaV_7 NsigmaV_9]);
% NsigmaB = min([NsigmaB_6 NsigmaB_7 NsigmaB_9]);
NsigmaV = min(cat(2, NsigmaV_6, NsigmaV_7, NsigmaV_9), [], 2);
NsigmaB = min(cat(2, NsigmaB_6, NsigmaB_7, NsigmaB_9), [], 2);
end

function D_nb26 = GetNeighbors(D, inPLX)
% 找D的26邻域，存在N×27表里，N为D的像素数量
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
% 只计算inPLX指示的部分
D_nb26 = D_neighbors(inPLX, :);
end

function [NsigmaV, NsigmaB] = GetNsigma2020(costheta, dI, NBS, D_nb26_foreground, D_nb26_background, epsilon1, epsilon2)
% 2020 lina版本的Nsigma求法，即\scriptP_C_(m,j)(x_s, x_s')
% INPUT
% costheta：N×27矩阵，每个体素与它26邻域体素的夹角的cos（PS: 本体素排在中间，第14个）
% Iout_nb26：N×27矩阵，每个体素与它26邻域体素的多通道滤波结果
% NBS：NBS旋转的数量×NBS像素数量矩阵，包括这个NBS类的所有旋转类型
% D_nb26_foreground：N×27矩阵，每个体素与它26邻域体素为前景标记为1，否则为0
% D_nb26_background：N×27矩阵，每个体素与它26邻域体素为背景标记为1，否则为0
% epsilon1, epsilon2：和为1的正常数，PWS-potential两个元素|▽I| / (|▽Imax|)和Similarity的权

[NBS_num, ~] = size(NBS); % 这个class NBS的数量
[lenPLX, ~] = size(costheta); % 所有像素数量
% NsigmaV = zeros(lenPLX, 1);
% NsigmaB = zeros(lenPLX, 1);
nsigmaV = zeros(lenPLX, NBS_num);
nsigmaB = zeros(lenPLX, NBS_num);
for j = 1:NBS_num
    costheta_NBSj = costheta(:,NBS(j,:));  % 取形状NBS的部分计算
    S = 1- abs(costheta_NBSj); % similarity，两个角度夹角↓→costheta↑→S↓
    e1 = exp(S); % 我发誓论文里没这步，而且【这不科学，这样S的范围就是1-e，而不是0-1，直接导致rho的值肯定大于1，那么1-rho就小于0了】
    
    dI_NBSj = dI(:,NBS(j,:)); % |▽I| ，强度差↓，dV↓
    dV_norm = dI_NBSj ./ (max(dI_NBSj, [], 2) + eps);  % |▽I| / (|▽Imax|)
    % 原来↓
    rho = epsilon1 * dV_norm + epsilon2 * e1; % 夹角↓，强度差↓，rho↓
    
%     NsigmaV = min(nsigmaV, sum(D_nb26_background(:,NBS(j,:)) .* rho, 2)); % xs, xs'同为L_V→0，xs, xs'不同类→rho
%     NsigmaB = min(nsigmaB, sum(D_nb26_foreground(:,NBS(j,:)) .* rho + D_nb26_background(:,NBS(j,:)) .* abs(1 - rho), 2)); %【1 - rho的范围根本不是0-1，这怎么能这样！】
    nsigmaV(:,j) = sum(D_nb26_background(:,NBS(j,:)) .* rho, 2); % xs, xs'同为L_V→0，xs, xs'不同类→rho
    nsigmaB(:,j) = sum(D_nb26_foreground(:,NBS(j,:)) .* rho + D_nb26_background(:,NBS(j,:)) .* abs(1 - rho), 2); %【1 - rho的范围根本不是0-1，这怎么能这样！】
    % 原来↑
%     % 新↓这里交换了rho中similarity和|▽I|的权重！
%     rho = 0 * dV_norm + 1 * e1; % 夹角↓，强度差↓，rho↓
%     nsigmaV(:,j) = sum(D_nb26_background(:,NBS(j,:)) .* rho, 2); % xs, xs'同为L_V→0，xs, xs'不同类→rho
%     rho = 0.5 * dV_norm + 0.5 * e1; % 夹角↓，强度差↓，rho↓
%     nsigmaB(:,j) = sum(D_nb26_foreground(:,NBS(j,:)) .* rho + D_nb26_background(:,NBS(j,:)) .* abs(1 - rho), 2); %【1 - rho的范围根本不是0-1，这怎么能这样！】
%     % 新↑
end
NsigmaV = min(nsigmaV, [], 2);
NsigmaB = min(nsigmaB, [], 2);
end

function [NsigmaV, NsigmaB] = GetNsigma2021(pVy, pBy, NBS, D_nb26_foreground, D_nb26_background)
% 2020 lina版本的Nsigma求法，即\scriptP_C_(m,j)(x_s, x_s')
% INPUT
% pVy：N×27矩阵，每个体素与它26邻域体素的pVy（PS: 本体素排在中间，第14个）
% pBy：N×27矩阵，每个体素与它26邻域体素的pBy
% NBS：NBS旋转的数量×NBS像素数量矩阵，包括这个NBS类的所有旋转类型
% D_nb26_foreground：N×27矩阵，每个体素与它26邻域体素为前景标记为1，否则为0
% D_nb26_background：N×27矩阵，每个体素与它26邻域体素为背景标记为1，否则为0
% epsilon1, epsilon2：和为1的正常数，PWS-potential两个元素|▽I| / (|▽Imax|)和Similarity的权

[NBS_num, ~] = size(NBS); % 这个class NBS的数量
[lenPLX, ~] = size(pVy); % 所有像素数量
% NsigmaV = zeros(lenPLX, 1);
% NsigmaB = zeros(lenPLX, 1);
nsigmaV = zeros(lenPLX, NBS_num);
nsigmaB = zeros(lenPLX, NBS_num);
for j = 1:NBS_num
    pV = pVy(:,NBS(j,:)) ./ (pVy(:,NBS(j,:)) + pBy(:,NBS(j,:)));
    % 原来↓
    rho = pV; % pV↓，rho↓
   nsigmaV(:,j) = sum(D_nb26_background(:,NBS(j,:)) .* rho, 2); % xs, xs'同为L_V→0，xs, xs'不同类→rho
    nsigmaB(:,j) = sum(D_nb26_foreground(:,NBS(j,:)) .* rho + D_nb26_background(:,NBS(j,:)) .* abs(1 - rho), 2); %【1 - rho的范围根本不是0-1，这怎么能这样！】
end
NsigmaV = min(nsigmaV, [], 2);
NsigmaB = min(nsigmaB, [], 2);
end