# script that computes for a given PPI network and functional PSN significances of overlaps

import sys, os
import numpy
import pandas

if __name__ == '__main__':

	network_name = sys.argv[1]
	fctPSN_file = sys.argv[2]
	fctPSN_name = sys.argv[3]
	num_rand = int(sys.argv[4])
	# either RAI or Jaccard_similarity
	edge_weight = sys.argv[5]
	min_pcc = float(sys.argv[6])
	max_pcc = float(sys.argv[7])

	path = '../data/katjas_data/PSN_analysis/'
	network_path = path

	# read in real PPI network
	if '_PSN' in network_name:
		nw_df = pandas.read_table(network_path + network_name + '.txt',header=0,usecols=['gene_a','gene_b',edge_weight])
	else:
		nw_df = pandas.read_table(network_path + network_name + '.txt',header=0,usecols=['gene_a','gene_b'])
		nw_df[edge_weight] = pandas.Series([1 for i in range(nw_df.shape[0])],index=nw_df.index)

	# generate lookup dictionary for gene IDs -> index in fctPSN matrix
	file1 = open(path + fctPSN_file + '.txt','r')
	for line in file1:
		fctPSN_gene_ids = line[:-1].split('\t')[1:]
		break
	file1.close()
	fctPSN_gene_id_dict = dict([(fctPSN_gene_ids[i],i) for i in range(len(fctPSN_gene_ids))])

	# read in fctPSN as matrix
	print('read in PCC matrix')
#	fctPSN_matrix = numpy.loadtxt(path + fctPSN_file,dtype='float',delimiter='\t',skiprows=1,usecols=range(1,len(fctPSN_gene_ids)+1))
	fctPSN_matrix = numpy.load(path + fctPSN_file + '.npy')

	# determine for every cutoff the number of possible protein pairs in the fctPSN_matrix between proteins that
	# are in the network
	# get all genes that are in the PPI network and in the fct PSN
	print('determine common pairs')
	overlap_genes = list(set(pandas.concat([nw_df['gene_a'],nw_df['gene_b']])).intersection(set(fctPSN_gene_ids)))
	# get the PCCs for all gene pair combinations and store in a list
	pccs_pos_pairs = []
	for i,gene_a in enumerate(overlap_genes):
		for gene_b in overlap_genes[i+1:]:
			pccs_pos_pairs.append(fctPSN_matrix[fctPSN_gene_id_dict[gene_a],fctPSN_gene_id_dict[gene_b]])

	# compute overlaps between real PPI network and fctPSN at various cutoffs
	num_cutoffs = 10
	step_size = (max_pcc - min_pcc)/float(num_cutoffs)
	cutoffs = [min_pcc + step_size*(i+1) for i in range(num_cutoffs)]
	# matrix that stores the sums of the edge weights for the edges in the overlap
	# first row is for real network
	# all other rows are for the randomized networks
	sum_matrix = numpy.zeros((num_rand+1,num_cutoffs),dtype=float)
	# add a column to the dataframe that lists for every PPI the PCC of the fct PSN if both genes are in the fct PSN
	nw_df['PCCs'] = nw_df.apply(lambda x: fctPSN_matrix[fctPSN_gene_id_dict[x['gene_a']],fctPSN_gene_id_dict[x['gene_b']]]
										  if (x['gene_a'] in fctPSN_gene_id_dict and x['gene_b'] in fctPSN_gene_id_dict)
										  else numpy.NaN,
								axis=1)
	last_cutoff = cutoffs[len(cutoffs)-1]
	sum_weights = nw_df.loc[nw_df['PCCs'] >= last_cutoff, [edge_weight]].sum()[edge_weight]
	sum_matrix[0,len(cutoffs)-1] = sum_weights
	# go backwards through the cutoffs, find all PCCs between the current and last cutoff and sum those up
	# add that sum to the sum from the previous cutoff and save in the matrix
	for i,cutoff in enumerate(cutoffs[:-1][::-1]):
		# select all rows with a PCC between the actual cutoff and the last cutoff
		sum_weights += nw_df.loc[(nw_df['PCCs'] >= cutoff) & (nw_df['PCCs'] < last_cutoff), [edge_weight]].sum()[edge_weight]
		sum_matrix[0,len(cutoffs)-i-2] = sum_weights
		last_cutoff = cutoff

	# compute overlaps for every random network
	for r in range(num_rand):

		print(r)

		# read in the network
		if '_PSN' in network_name:
			rn_df = pandas.read_table(network_path + network_name + '/' + network_name + '_rand_network_' + str(r) + '.txt',header=0,usecols=['gene_a','gene_b',edge_weight])
		else:
			rn_df = pandas.read_table(network_path + network_name + '/' + network_name + '_rand_network_' + str(r) + '.txt',header=0,usecols=['gene_a','gene_b'])
			rn_df[edge_weight] = pandas.Series([1 for i in range(rn_df.shape[0])],index=rn_df.index)

		# add a column to the dataframe that lists for every PPI the PCC of the fct PSN if both genes are in the fct PSN
		rn_df['PCCs'] = rn_df.apply(lambda x: fctPSN_matrix[fctPSN_gene_id_dict[x['gene_a']],fctPSN_gene_id_dict[x['gene_b']]]
											  if (x['gene_a'] in fctPSN_gene_id_dict and x['gene_b'] in fctPSN_gene_id_dict)
											  else False,
									axis=1)
		last_cutoff = cutoffs[len(cutoffs)-1]
		sum_weights = rn_df.loc[rn_df['PCCs'] >= last_cutoff, [edge_weight]].sum()[edge_weight]
		sum_matrix[r+1,len(cutoffs)-1] = sum_weights
		# go backwards through the cutoffs, find all PCCs between the current and last cutoff and sum those up
		# add that sum to the sum from the previous cutoff and save in the matrix
		for i,cutoff in enumerate(cutoffs[:-1][::-1]):
			# select all rows with a PCC between the actual cutoff and the last cutoff
			sum_weights += rn_df.loc[(rn_df['PCCs'] >= cutoff) & (rn_df['PCCs'] < last_cutoff), [edge_weight]].sum()[edge_weight]
			sum_matrix[r+1,len(cutoffs)-i-2] = sum_weights
			last_cutoff = cutoff

	# write out results
	outfile = path + 'nw_overlaps_per_PCCcutoffs_sum_edge_weights_' + edge_weight + '_' + network_name + '_' + fctPSN_name + '.txt'
	target = open(outfile,'w')
	target.write('cutoff\tsum_edge_weights_real\tmean_sum_edge_weights_rand\tstd_sum_edge_weights_rand\tmedian_sum_edge_weights_rand\t' + \
				 '2.5perc_sum_edge_weights_rand\t97.5perc_sum_edge_weights_rand\tz_score\tpvalue\tnum_rand_ge\tnum_pairs_common_space\n')
	for i,cutoff in enumerate(cutoffs):
		real_sum = sum_matrix[0,i]
		rand_mean = numpy.mean(sum_matrix[1:,i])
		rand_std = numpy.std(sum_matrix[1:,i])
		if rand_std == 0:
			zscore = numpy.NaN
		else:
			zscore = (real_sum - rand_mean)/rand_std
		num_rand_ge = len(filter(lambda v: v >= sum_matrix[0,i],sum_matrix[1:,i]))
		pvalue = num_rand_ge/float(num_rand)
		num_pos_pairs = len(filter(lambda v: v >= cutoff,pccs_pos_pairs))
		target.write(str(cutoff) + '\t' + str(real_sum) + '\t' + str(rand_mean) + '\t' + \
					 str(rand_std) + '\t' + str(numpy.median(sum_matrix[1:,i])) + '\t' + \
					 str(numpy.percentile(sum_matrix[1:,i],2.5)) + '\t' + str(numpy.percentile(sum_matrix[1:,i],97.5)) + '\t' + \
					 str(zscore) + '\t' + str(pvalue) + '\t' + str(num_rand_ge) + '\t' + str(num_pos_pairs) + '\n')
	target.close()

	outfile = path + 'raw_nw_overlaps_per_PCCcutoffs_sum_edge_weights_' + edge_weight + '_' + network_name + '_' + fctPSN_name + '.npy'
	numpy.save(outfile,sum_matrix)
