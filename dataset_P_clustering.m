%{
Author: S. Watanabe
Last updated 12/28/2024

This program is for the clustering gene expression profiles in ASD

Input data are from Parikshak et al. Nature (2016)

variables
p_symbol	gene symbol
p_exp		gene expression
p_ctl_i		control sample index (age range 10-60 years)
p_asd_i		ASD sample index (age range 10-60 years)
%}


%list of cell type-specific genes from McKenzie SciRep (2018)
ct_filename = '41598_2018_27293_MOESM2_ESM.xlsx';
ct_sheetname = 'top_human_enrich';
ct_gene_tbl = readtable(ct_filename, Sheet=ct_sheetname, NumHeaderLines=2);

celltype_list = {'ast', 'end', 'mic', 'neu', 'oli'};
n_celltypes = length(celltype_list);

ct_gene = {};
for ct=1:n_celltypes
	ct_range = (ct-1)*1000 + (1:1000);
	ct_gene{ct} = ct_gene_tbl.gene(ct_range);
end


%normalized expression in dataset P
p_exp_norm = (p_exp - mean(p_exp(:, p_ctl_i), 2)) ./ std(p_exp(:, p_ctl_i), 0, 2);

%number of cell type specific genes used for clustering
n_ct_genes_clust = 50;

%number of clusters
n_asd_grp = 3;


p_ct_gene = cell(1,5);
p_ct_i = cell(1,5);
p_ct_i_cl = cell(1,5);

for ct=1:5
	[~, i_all] = ismember(ct_gene{ct}, p_symbol);
	p_ct_i{ct} = i_all(i_all>0);
	p_ct_gene{ct} = p_symbol(p_ct_i{ct});
	p_ct_i_cl{ct} = p_ct_i{ct}(1:n_ct_genes_clust);
end
p_ct_i_cl_4ct = [p_ct_i_cl{1}; p_ct_i_cl{3}; p_ct_i_cl{4}; p_ct_i_cl{5}];


%clustering
clust_data_clust = p_exp_norm(p_ct_i_cl_4ct, p_asd_i);
tree_s = linkage(clust_data_clust', 'ward'); 
[~, ~, perm_s] = dendrogram(tree_s, 0);
p_grp_raw = cluster_tree(tree_s, n_asd_grp);


%group with highest mic expression = 1, lowest neu expression = 2
grp_mic_mean = zeros(1, n_asd_grp);
grp_neu_mean = zeros(1, n_asd_grp);
for i=1:n_asd_grp
	grp_mic_mean(i) = mean(mean(p_exp(p_ct_i_cl{3}, p_asd_i(p_grp_raw==i))));
	grp_neu_mean(i) = mean(mean(p_exp(p_ct_i_cl{4}, p_asd_i(p_grp_raw==i))));
end
[~, grp_mic_mean_max] = max(grp_mic_mean);
[~, grp_neu_mean_max] = min(grp_neu_mean);

p_grp = 3 * ones(1, length(p_asd_i));
p_grp(p_grp_raw==grp_mic_mean_max) = 1;
p_grp(p_grp_raw==grp_neu_mean_max) = 2;

