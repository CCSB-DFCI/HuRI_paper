# script that computes for a given PPI network and functional PSN significances of overlaps

import sys, os
import numpy
import pandas


# read in all the co annotations into a list of dicts
# read in the network, determine whether it is a PPI or PSN to get the edge weights right
# for the real network and every co-annotation dataset, determine the sum of edge weights of all edges that share an annotation
# and keep track of common space: how many proteins in network have an annotation
# save in list of lists the edge weight sums for all the random networks

if __name__ == '__main__':

	network_name = sys.argv[1]
	num_rand = int(sys.argv[2])
	# either RAI or Jaccard_similarity
	edge_weight = sys.argv[3]

	path = '../data/katjas_data/PSN_analysis/'
	network_path = path
	annot_path = path
	out_path = path

	annot_dataset_names = ['Bioplex2.0_complex_info','cell_atlas_localization_info','reactome_pathway_info']
	annot_datasets = []
	for name in annot_dataset_names:
		file1 = open(annot_path + name + '.txt','r')
		entries = file1.readlines()
		entries = [line[:-1].split('\t') for line in entries[1:]]
		file1.close()
		annot_dict = dict([(t[0],set(t[1].split('|'))) for t in entries])
		annot_datasets.append(annot_dict)

	# read in real PPI network
	if '_PSN' in network_name:
		nw_df = pandas.read_table(network_path + network_name + '.txt',header=0,usecols=['gene_a','gene_b',edge_weight])
	else:
		nw_df = pandas.read_table(network_path + network_name + '.txt',header=0,usecols=['gene_a','gene_b'])
		nw_df[edge_weight] = pandas.Series([1 for i in range(nw_df.shape[0])],index=nw_df.index)

	# add to network dataframe a column for each annot dataset with 1 when both genes share annotation, 0 if they don't and nan if
	# at least one of the two genes has no annotation available
	# calculate the sum of the edge weights of all edges that have a co-annotation
	sum_matrix = numpy.zeros((num_rand+1,len(annot_dataset_names)),dtype=float)
	for d,d_name in enumerate(annot_dataset_names):
		nw_df[d_name] = nw_df.apply(lambda x: 0
											  if x['gene_a'] in annot_datasets[d] and x['gene_b'] in annot_datasets[d]
											  else numpy.nan,
									axis=1)
		nw_df[d_name] = nw_df.apply(lambda x: 1
											  if x[d_name] == 0 and len(annot_datasets[d][x['gene_a']].intersection(annot_datasets[d][x['gene_b']])) > 0
											  else x[d_name],
									axis=1)
		sum_matrix[0,d] = nw_df.loc[nw_df[d_name] == 1,[edge_weight]].sum()[edge_weight]


	# compute overlaps for every random network
	for r in range(num_rand):

		print(r)

		# read in the network
		if '_PSN' in network_name:
			rn_df = pandas.read_table(network_path + network_name + '/' + network_name + '_rand_network_' + str(r) + '.txt',header=0,usecols=['gene_a','gene_b',edge_weight])
		else:
			rn_df = pandas.read_table(network_path + network_name + '/' + network_name + '_rand_network_' + str(r) + '.txt',header=0,usecols=['gene_a','gene_b'])
			rn_df[edge_weight] = pandas.Series([1 for i in range(rn_df.shape[0])],index=rn_df.index)

		# add annotation info
		for d,d_name in enumerate(annot_dataset_names):
			rn_df[d_name] = rn_df.apply(lambda x: 0
												  if x['gene_a'] in annot_datasets[d] and x['gene_b'] in annot_datasets[d]
												  else numpy.nan,
										axis=1)
			rn_df[d_name] = rn_df.apply(lambda x: 1
												  if x[d_name] == 0 and len(annot_datasets[d][x['gene_a']].intersection(annot_datasets[d][x['gene_b']])) > 0
												  else x[d_name],
										axis=1)
			sum_matrix[r+1,d] = rn_df.loc[rn_df[d_name] == 1,[edge_weight]].sum()[edge_weight]


	# write out results
	outfile = out_path + 'nw_overlaps_coAnnot_sum_edge_weights_' + edge_weight + '_' + network_name + '.txt'
	target = open(outfile,'w')
	target.write('fct_dataset\tsum_edge_weights_real\tmean_sum_edge_weights_rand\tstd_sum_edge_weights_rand\tmedian_sum_edge_weights_rand\t' + \
				 '2.5perc_sum_edge_weights_rand\t97.5perc_sum_edge_weights_rand\tz_score\tpvalue\tnum_rand_ge\tnum_pairs_common_space\n')
	for d,d_name in enumerate(annot_dataset_names):
		real_sum = sum_matrix[0,d]
		rand_mean = numpy.mean(sum_matrix[1:,d])
		rand_std = numpy.std(sum_matrix[1:,d])
		if rand_std == 0:
			zscore = numpy.NaN
		else:
			zscore = (real_sum - rand_mean)/rand_std
		num_rand_ge = len(filter(lambda v: v >= sum_matrix[0,d],sum_matrix[1:,d]))
		pvalue = num_rand_ge/float(num_rand)
		num_pos_pairs = nw_df.loc[(nw_df[d_name] == 0) | (nw_df[d_name] == 1),].shape[0]
		target.write(d_name + '\t' + str(real_sum) + '\t' + str(rand_mean) + '\t' + \
					 str(rand_std) + '\t' + str(numpy.median(sum_matrix[1:,d])) + '\t' + \
					 str(numpy.percentile(sum_matrix[1:,d],2.5)) + '\t' + str(numpy.percentile(sum_matrix[1:,d],97.5)) + '\t' + \
					 str(zscore) + '\t' + str(pvalue) + '\t' + str(num_rand_ge) + '\t' + str(num_pos_pairs) + '\n')
	target.close()
