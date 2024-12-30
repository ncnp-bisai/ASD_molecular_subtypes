%{
Author: S. Watanabe
Last updated 12/28/2024

This program is for the clustering gene expression profiles in ASD

Input data are from Zhang et al. PNAS (2023)

variables
z_symbol	gene symbol
z_exp		gene expression
z_asd_i		ASD sample index (age range 10-60 years)
z_ctl_i		control sample index (age range 10-60 years)
%}


%normalized expression in dataset Z
z_exp_norm = (z_exp - mean(z_exp(:, z_ctl_i), 2)) ./ std(z_exp(:, z_ctl_i), 0, 2);

%number of cell type specific genes used for clustering
n_ct_genes_clust = 50;


[p_z_symbol, p_i, z_i] = intersect(p_symbol, z_symbol);
p_z_exp_norm = p_exp_norm(p_i, :);
z_p_exp_norm = z_exp_norm(z_i, :);

p_z_ct_gene = cell(1,5);
p_z_ct_i = cell(1,5);
p_z_ct_i_cl = cell(1,5);

for ct=1:5
	[~, i_all] = ismember(ct_gene{ct}, p_z_symbol);
	p_z_ct_i{ct} = i_all(i_all>0);
	p_z_ct_gene{ct} = p_z_symbol(p_z_ct_i{ct});
	p_z_ct_i_cl{ct} = p_z_ct_i{ct}(1:n_ct_genes_clust);
end
p_z_ct_i_cl_4ct = [p_z_ct_i_cl{1}; p_z_ct_i_cl{3}; p_z_ct_i_cl{4}; p_z_ct_i_cl{5}];


%h-sch correlation
p_z_corr_data = p_z_exp_norm(p_z_ct_i_cl_4ct, p_asd_i) - mean(p_z_exp_norm(p_z_ct_i_cl_4ct, p_asd_i), 2);
z_p_corr_data = z_p_exp_norm(p_z_ct_i_cl_4ct, z_asd_i) - mean(z_p_exp_norm(p_z_ct_i_cl_4ct, z_asd_i), 2);
p_z_corr = corr(p_z_corr_data, z_p_corr_data);


%grouping based on average p-z correlation 
p_z_corr_grpav = [mean(p_z_corr(p_asd_grp==1, :), 1); ...
					mean(p_z_corr(p_asd_grp==2, :), 1); ...
					mean(p_z_corr(p_asd_grp==3, :), 1)];
[~, z_grp] = max(p_z_corr_grpav, [], 1);

