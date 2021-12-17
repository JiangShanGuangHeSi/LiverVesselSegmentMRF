function [Idx,Ctrs] = Kmean_reorder(idx,ctrs)
% �������壺��idx��ctrs�е����ݰ���ctrs��С�����˳����������
K = length(ctrs);
Ctrs_index = [ctrs (1:length(ctrs))'];
sort_Ctrs = sortrows(Ctrs_index,1);% sort_Ctrs(:,1)�ɵ����ߴ�ž������ģ�sort_Ctrs(:,2)��Ŷ�Ӧ��ԭʼ����k
K_index = cell(1,K);
for k = 1:K % ��idx��ԭʼ�����k�ֱ�洢��K_index�ṹ��
    K_index{k} = find(idx==k);
end
Idx = zeros(size(idx));
Ctrs = zeros(size(ctrs));
for i = 1:K
    Ctrs(i) = sort_Ctrs(i,1);
    Idx(K_index{sort_Ctrs(i,2)}) = i;
end
end