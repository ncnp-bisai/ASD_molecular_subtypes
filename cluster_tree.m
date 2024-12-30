function data_cluster = cluster_tree(tree, n_clusters)

n_br = size(tree, 1);
n_data = n_br +1;

tree_e = [(n_data+1:2*n_data-1)', tree];

max_dist = tree_e(n_br-n_clusters+1, 4);

top_node = zeros(1, n_data);

for i=1:n_data
	node = i;
	
	while 1
		row = find(node==tree_e(:,2) | node==tree_e(:,3));
		if tree_e(row, 4) > max_dist
			break
		end
		node = tree_e(row, 1);
	end
	
	top_node(i) = node;
end

top_node_list = unique(top_node);

data_cluster = zeros(1, n_data);
for i=1:n_data
	data_cluster(i) = find(top_node(i)==top_node_list);
end

