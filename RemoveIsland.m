function Iout = RemoveIsland(Iin, min_size)
% 输入二值图像Iin，要保留的最小连通域大小min_size
% 输出去除小于min_size连通域的二值图像

I_binery = Iin ~= 0;
sIdx = regionprops3(I_binery, 'VoxelIdxList'); % sIdx类型为tabel

Iout = zeros(size(Iin));

for i = 1:height(sIdx) % 未定义与 'table' 类型的输入参数相对应的函数 'length'。请改用 HEIGHT、WIDTH 或 SIZE 函数。
%     num_sIdx = length(sIdx(i, 1).VoxelIdxList);
    island_idx = sIdx.VoxelIdxList(i); % 这是一个cell
    island_idx = island_idx{1,1}; % cell中的第1,1个东西，一个表
    num_sIdx = length(island_idx);
    if(num_sIdx > min_size)
        Iout(island_idx) = 1;
    end
end

end