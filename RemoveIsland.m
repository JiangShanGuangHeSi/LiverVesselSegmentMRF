function Iout = RemoveIsland(Iin, min_size)
% �����ֵͼ��Iin��Ҫ��������С��ͨ���Сmin_size
% ���ȥ��С��min_size��ͨ��Ķ�ֵͼ��

I_binery = Iin ~= 0;
sIdx = regionprops3(I_binery, 'VoxelIdxList'); % sIdx����Ϊtabel

Iout = zeros(size(Iin));

for i = 1:height(sIdx) % δ������ 'table' ���͵�����������Ӧ�ĺ��� 'length'������� HEIGHT��WIDTH �� SIZE ������
%     num_sIdx = length(sIdx(i, 1).VoxelIdxList);
    island_idx = sIdx.VoxelIdxList(i); % ����һ��cell
    island_idx = island_idx{1,1}; % cell�еĵ�1,1��������һ����
    num_sIdx = length(island_idx);
    if(num_sIdx > min_size)
        Iout(island_idx) = 1;
    end
end

end