# script to test whether TiP genes are close in the networks
# script will test different ways of determining TiP genes, different networks and
# different ways of determining closeness in networks
# it will compute these tests for every tissue separately

import sys, os, pickle, random
import numpy, pandas, igraph
import database_utils, random_graphs, config


# function that returns a graph object of the largest connected component
# formed by the given set of seed genes and graph
def get_LCC(graph,seed_genes):

	sg = graph.subgraph(seed_genes)
	clusters = sg.clusters(mode='WEAK')
	return clusters.giant()


# function that computes for a given set of seed genes the average shortest path length
# and returns it
def get_avg_shortest_paths(graph,seed_genes):

	paths = []
	shortest_paths_lists = graph.shortest_paths_dijkstra(source=seed_genes,target=seed_genes)
	for shortest_paths in shortest_paths_lists:
		shortest_paths.remove(0)
		if len(shortest_paths) > 0:
			min_path = min(shortest_paths)
			paths.append(min_path)
	if len(paths) > 0:
		asp = (numpy.sum(paths))/float(len(paths))
	else:
		asp = 0
	return asp


# function that returns the number of edges linking seed genes directly to each other in
# the given graph
def get_num_edges(graph,seed_genes):

	sg = graph.subgraph(seed_genes)
	return len(sg.es)


# function that returns the non-redundant count of interaction partners for the
# given gene set
def get_num_interaction_partners(graph,gene_set):

	neighbor_lists = graph.neighborhood(gene_set)
	# neighborhood returns in every list of neighbors for a given node the node itself at the first
	# position of the neighbor_list -> to be removed
	all_neighbors = [n for n_list in neighbor_lists for n in n_list[1:]]
	all_neighbors = set([graph.vs[n]['name'] for n in all_neighbors])

	return all_neighbors


# returns a CSN for the given source network and tissue
def get_CSN_for_tissue_and_source_network(name_source_network,tissue):

	file1 = open(config.analysis_path + name_source_network + '/GTExv6_' + name_source_network + '_CSNs.txt','r')
	entries = file1.readlines()
	file1.close()
	tissues = entries[0][:-1].split('\t')
	edge_list = []
	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		gene_a = tab_list[0]
		gene_b = tab_list[1]
		frac_samples = float(tab_list[tissues.index(tissue)])
		if frac_samples > config.edge_expr_cutoff:
			edge_list.append((gene_a,gene_b))
	nodes = set()
	for edge in edge_list:
		nodes.add(edge[0])
		nodes.add(edge[1])
	CSN = igraph.Graph()
	CSN.add_vertices(list(nodes))
	CSN.add_edges(edge_list)

	return CSN


def test_TiP_network_closeness(TiPmatrix,TiP_cutoff,closeness_types,CSNs,network_name,num_rand):

	conditions = list(CSNs.keys())
	conditions.sort()

	# create output files
	targets = []
	for closeness_type in closeness_types:
		if 'testis' not in conditions:
			outfile = config.analysis_path + network_name + '/GTExv6_NT_' + network_name + '_TiPgenes_' + str(TiP_cutoff) + '_' + closeness_type + '_closeness_test_summary.txt'
		else:
			outfile = config.analysis_path + network_name + '/GTExv6_' + network_name + '_TiPgenes_' + str(TiP_cutoff) + '_' + closeness_type + '_closeness_test_summary.txt'
		target = open(outfile,'w')
		target.write('condition\treal_obs\tnum_rand\tnum_rand_obs_greater\tnum_rand_obs_smaller\t' + \
					'mean_rand_obs\tstd_rand_obs\tmedian_rand_obs\t5th_perc_rand_obs\t95_perc_rand_obs\t' + \
					'z_score\tp_value_greater\tp_value_smaller\tnum_TiP_genes_in_CSN\tnum_PPI_TiP_genes_in_CSN\n')
		targets.append(target)

	matrices = []
	for closeness_type in closeness_types:
		matrices.append(numpy.empty(shape=(num_rand,len(conditions)),dtype=float))

	# for every tissue
	for c,condition in enumerate(conditions):

		# determine the TiP genes for this tissue
		print(condition)
		CSN = CSNs[condition]
		TiP_genes = set(TiPmatrix.loc[TiPmatrix[condition]>=TiP_cutoff,].index.tolist())
		CSN_nodes = set(CSN.vs['name'])
		# get a graph that is the LCC of the CSN in which I can only determine avg shortest path of TiP genes
		LCC_CSN = CSN.clusters(mode='WEAK').giant()
		LCC_CSN_nodes = set(LCC_CSN.vs['name'])
		TiP_genes_in_CSN = TiP_genes.intersection(CSN_nodes)
		TiP_genes_in_LCC_CSN = TiP_genes.intersection(LCC_CSN_nodes)
		real_obs_list = []
		# determine for every closeness type the real observation and save in list
		for closeness_type in closeness_types:
			if closeness_type == 'nl':
				if len(TiP_genes_in_CSN) < 2:
					real_obs = 'NaN'
				else:
					real_obs = get_num_edges(CSN,TiP_genes_in_CSN)
			elif closeness_type == 'lcc':
				if len(TiP_genes_in_LCC_CSN) < 2:
					real_obs = 'NaN'
				else:
					lcc_tip = get_LCC(LCC_CSN,TiP_genes_in_LCC_CSN)
					real_obs = len(lcc_tip.vs)
			elif closeness_type == 'asp':
				if len(TiP_genes_in_LCC_CSN) < 2:
					real_obs = 'NaN'
				else:
					real_obs = get_avg_shortest_paths(LCC_CSN,TiP_genes_in_LCC_CSN)
			else:
				print('unknown type of closeness measure')
				sys.exit()
			real_obs_list.append(real_obs)

		print(closeness_types)
		print(real_obs_list)

		rand_obs_lists = [[] for i in range(0,len(closeness_types))]
		# determine the random observations in all random networks for every closeness type and save them in lists
		if len(TiP_genes_in_LCC_CSN) > 1:
			for r in range(num_rand):
				rand_file = config.rand_dir + network_name + '/random_CSN_networks/' + condition + '/' + network_name + '_rand_network_' + str(r) + '.txt'
				rg = random_graphs.get_real_network(rand_file,None)
				TiP_in_rg = TiP_genes.intersection(set(rg.vs['name']))
				LCC_rg = rg.clusters(mode='WEAK').giant()
				TiP_in_LCC_rg = TiP_genes.intersection(set(LCC_rg.vs['name']))
				for t,closeness_type in enumerate(closeness_types):
					if closeness_type == 'nl':
						if len(TiP_in_rg) > 1:
							rand_nl = get_num_edges(rg,TiP_in_rg)
							rand_obs_lists[t].append(rand_nl)
						else:
							rand_obs_lists[t].append('NaN')
					elif closeness_type == 'lcc':
						if len(TiP_in_LCC_rg) > 1:
							lcc_rg_tip = get_LCC(LCC_rg,TiP_in_LCC_rg)
							rand_obs_lists[t].append(len(lcc_rg_tip.vs))
						else:
							rand_obs_lists[t].append('NaN')
					elif closeness_type == 'asp':
						if len(TiP_in_LCC_rg) > 1:
							rand_asp = get_avg_shortest_paths(LCC_rg,TiP_in_LCC_rg)
							rand_obs_lists[t].append(rand_asp)
						else:
							rand_obs_lists[t].append('NaN')

		for t,closeness_type in enumerate(closeness_types):
			matrices[t][:,c] = rand_obs_lists[t]

		print(rand_obs_lists)
		# calculate statistics and write out the results
		for t,closeness_type in enumerate(closeness_types):
			num_valid = len(list(filter(lambda v: str(v).isdigit(), rand_obs_lists[t])))
			# if there weren't enough TiPs for this tissue or in enough of the random networks,
			# write out this test as invalid
			if real_obs_list[t] == 'NaN' or num_valid < 100:
				targets[t].write(condition + '\tNaN\t' + str(num_rand) + '\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\t' + \
								 (str(len(TiP_genes_in_LCC_CSN)) if condition in ['asp','lcc'] else str(len(TiP_genes_in_CSN))) + '\t' + \
								 (str(len(get_num_interaction_partners(LCC_CSN,TiP_genes_in_LCC_CSN))) if \
								 condition in ['asp','lcc'] else str(len(get_num_interaction_partners(CSN,TiP_genes_in_CSN)))) + '\n')
			else:
				greater_rand_obs = list(filter(lambda l: l >= real_obs_list[t],rand_obs_lists[t]))
				smaller_rand_obs = list(filter(lambda l: l <= real_obs_list[t],rand_obs_lists[t]))
				p_value_g = len(greater_rand_obs)/float(num_valid)
				p_value_s = len(smaller_rand_obs)/float(num_valid)
				mean_rand_obs = numpy.mean(rand_obs_lists[t])
				std_rand_obs = numpy.std(rand_obs_lists[t])
				median_rand_obs = numpy.median(rand_obs_lists[t])
				perc5_rand_obs = numpy.percentile(rand_obs_lists[t],5)
				perc95_rand_obs = numpy.percentile(rand_obs_lists[t],95)
				if std_rand_obs == 0:
					z_score = 'NaN'
				else:
					z_score = (real_obs_list[t] - mean_rand_obs)/std_rand_obs
				targets[t].write(condition + '\t' + str(real_obs_list[t]) + '\t' + str(num_valid) + '\t' + str(len(greater_rand_obs)) + '\t' + \
							str(len(smaller_rand_obs)) + '\t' + str(mean_rand_obs) + '\t' + str(std_rand_obs) + '\t' + \
							str(median_rand_obs) + '\t' + str(perc5_rand_obs) + '\t' + str(perc95_rand_obs) + '\t' + \
							str(z_score) + '\t' + str(p_value_g) + '\t' + str(p_value_s) + '\t' + \
							 (str(len(TiP_genes_in_LCC_CSN)) if condition in ['asp','lcc'] else str(len(TiP_genes_in_CSN))) + '\t' + \
							 (str(len(get_num_interaction_partners(LCC_CSN,TiP_genes_in_LCC_CSN))) if condition in ['asp','lcc'] else str(len(get_num_interaction_partners(CSN,TiP_genes_in_CSN)))) + '\n')

	for target in targets:
		target.close()

	# write out random data
	for t,closeness_type in enumerate(closeness_types):
		if 'testis' not in conditions:
			outfile = config.analysis_path + network_name + '/GTExv6_NT_' + network_name + '_TiPgenes_' + str(TiP_cutoff) + '_' + closeness_type + '_rand_results.txt'
		else:
			outfile = config.analysis_path + network_name + '/GTExv6_' + network_name + '_TiPgenes_' + str(TiP_cutoff) + '_' + closeness_type + '_rand_results.txt'
		target = open(outfile,'w')
		target.write('\t'.join(conditions) + '\n')
		for i in range(num_rand):
			target.write('\t'.join([str(v) for v in matrices[t][i,]]) + '\n')
		target.close()


if __name__ == '__main__':


	if len(sys.argv) < 5:
		print(\
		"""insufficient number of arguments:
		network name
		closeness types comma separated, options: lcc, asp, nl
		TiP cutoff
		TiP gene determination with (True) or without (False) testis
		number of randomizations""")
		sys.exit()
	else:
		network_name = sys.argv[1]
		closeness_types = sys.argv[2].split(',')
		TiP_cutoff = float(sys.argv[3])
		with_testis = sys.argv[4]
		num_rand = int(sys.argv[5])

	print(network_name, TiP_cutoff, closeness_types, with_testis)

	print('load gene expression data')
	if with_testis == 'True':
		expr_file = config.expr_file
		pref_expr_file = config.pref_expr_file
	else:
		expr_file = config.NT_expr_file
		pref_expr_file = config.NT_pref_expr_file
	GTEx = pandas.read_table(expr_file,index_col=0)
	TiPmatrix = pandas.read_table(pref_expr_file,index_col=0)
	GTEx_tissues = TiPmatrix.columns
	for tissue in GTEx_tissues:
		TiPmatrix.loc[TiPmatrix.index.isin(GTEx.index[GTEx[tissue]<=config.gene_expr_cutoff]),[tissue]] = numpy.NaN

	### load the CSNs
	print('load CSNs')
	CSNs = {}
	for condition in GTEx_tissues:
		CSNs[condition] = get_CSN_for_tissue_and_source_network(network_name,condition)

	test_TiP_network_closeness(TiPmatrix,TiP_cutoff,closeness_types,CSNs,network_name,num_rand)
