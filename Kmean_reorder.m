function [Idx,Ctrs] = Kmean_reorder(idx,ctrs)
% 函数意义：将idx和ctrs中的内容按照ctrs由小到大的顺序重新排序
K = length(ctrs);
Ctrs_index = [ctrs (1:length(ctrs))'];
sort_Ctrs = sortrows(Ctrs_index,1);% sort_Ctrs(:,1)由低至高存放聚类中心；sort_Ctrs(:,2)存放对应的原始类别号k
K_index = cell(1,K);
for k = 1:K % 将idx中原始分类号k分别存储在K_index结构中
    K_index{k} = find(idx==k);
end
Idx = zeros(size(idx));
Ctrs = zeros(size(ctrs));
for i = 1:K
    Ctrs(i) = sort_Ctrs(i,1);
    Idx(K_index{sort_Ctrs(i,2)}) = i;
end
end