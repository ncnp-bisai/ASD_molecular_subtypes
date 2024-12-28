%{
Author: S. Watanabe
Last updated 12/28/2024

This program is for the clustering gene expression profiles in ASD

Input data are from Parikshak et al. Nature (2016)

variables
h_symbol	gene symbol
h_exp		gene expression
asd_i		ASD sample index (age range 10-60 years)
ctl_i		control sample index (age range 10-60 years)
%}

%normalized expression
h_exp_norm = (h_exp - mean(h_exp(:, ctl_i), 2)) ./ std(h_exp(:, ctl_i), 0, 2);

%number of cell type specific genes used for clustering
n_ct_genes_clust = 50;

%number of clusters
n_asd_grp = 3;

celltype_list = {'ast', 'end', 'mic', 'neu', 'oli'};
n_celltypes = length(celltype_list);


%cell type gene list from McKenzie SciRep (2018)
cg_gene_filename = '41598_2018_27293_MOESM2_ESM.xlsx';
sheetname = 'top_human_enrich';
ct_gene_tbl = readtable(cg_gene_filename, Sheet=sheetname, NumHeaderLines=2);

ct_gene = {};
for ct=1:n_celltypes
	ct_range = (ct-1)*1000 + (1:1000);
	ct_gene{ct} = ct_gene_tbl.gene(ct_range);
end


h_ct_gene = cell(1,5);
h_ct_i = cell(1,5);
h_ct_i_cl = cell(1,5);

for ct=1:5
	[~, i_all] = ismember(ct_gene{ct}, h_symbol);
	h_ct_i{ct} = i_all(i_all>0);
	h_ct_gene{ct} = h_symbol(h_ct_i{ct});
	h_ct_i_cl{ct} = h_ct_i{ct}(1:n_ct_genes_clust);
end
h_ct_i_cl_4ct = [h_ct_i_cl{1}; h_ct_i_cl{3}; h_ct_i_cl{4}; h_ct_i_cl{5}];

%clustering
clust_data_clust = h_exp_norm(h_ct_i_cl_4ct, asd_i);
tree_s = linkage(clust_data_clust', 'ward'); 
[~, ~, perm_s] = dendrogram(tree_s, 0);
asd_grp = cluster_tree(tree_s, n_asd_grp);


