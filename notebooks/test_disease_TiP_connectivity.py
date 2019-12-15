# script that assesses the significance of the causal genes for a given tissue to
# interact with TiP genes for that tissue for a given TiP cutoff

import sys, os
import numpy, pandas, igraph
import config, random_graphs, database_utils

def get_CG_tissue_pairs_for_TiP_cutoff_sel_diseases_expr_filtered(cursor,TiP_cutoff):

	query = """select distinct a.ensembl_gene_id,a.disease_tissue
			from (select *
				from hi2018_paper.OMIM_disease_tissue
				where number_tissues<=3 and gene_type='protein_coding' and tissue_expr_value>5 and
					disease_tissue!='testis' and tissue_TiP_value<{}) as a
			where a.omim_id not in
				(select distinct omim_id
				from hi2018_paper.OMIM_disease_tissue
				where number_tissues<=3 and gene_type='protein_coding' and tissue_expr_value>5 and
					tissue_TiP_value>={} and disease_tissue!='testis')""".format(TiP_cutoff,TiP_cutoff)

	cursor.execute(query)
	tissue_CG_dict = {}
	for row in cursor:
		gene_id = row[0]
		tissue = row[1]
		if tissue not in tissue_CG_dict:
			tissue_CG_dict[tissue] = set()
		tissue_CG_dict[tissue].add(gene_id)

	return tissue_CG_dict


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



if __name__ == '__main__':

	network_name = sys.argv[1]
	num_rand = int(sys.argv[2])
	TiP_cutoff = float(sys.argv[3])

	connect = database_utils.get_connection()
	cursor = connect.cursor()

	# get the causal gene hit tissue assignments
	print('load the disease gene - tissue assignments')
	tissue_CG_dict = get_CG_tissue_pairs_for_TiP_cutoff_sel_diseases_expr_filtered(cursor,TiP_cutoff)
	CG_tissues = list(tissue_CG_dict.keys())
	CG_tissues.sort()

	# read in the gene expression data
	print('load the GTEx gene expression data')
	GTEx = pandas.read_table(config.NT_expr_file,index_col=0)
	TiPmatrix = pandas.read_table(config.NT_pref_expr_file,index_col=0)
	GTEx_tissues = TiPmatrix.columns
	for tissue in GTEx_tissues:
		TiPmatrix.loc[TiPmatrix.index.isin(GTEx.index[GTEx[tissue]<=5]),[tissue]] = numpy.NaN

	# get the CSNs
	print('load the CSNs')
	CSNs = {}
	for condition in GTEx_tissues:
		CSNs[condition] = get_CSN_for_tissue_and_source_network(network_name,condition)

	outfile = '../data/katjas_data/tissue_spec_disease_edgotyping/' + network_name + '_CG_TiP_connect_CSNs_TiPcutoff' + str(TiP_cutoff) + '.txt'
	target = open(outfile,'w')
	target.write('tissue\tnum_TiP_genes_in_CSN\tnum_CGs_in_CSN\tnum_rand\t' + \
	 			'num_CG_TiP_PPIs\tnum_rand_CG_TiP_PPIs_as_extreme\tPPI_pvalue\tPPI_zscore\tPPI_5perc\tPPI_median\tPPI_95perc\tPPI_mean\tPPI_std\t' + \
				'num_CGs_with_TiP_PPIs\tnum_rand_CGs_with_TiP_PPIs_as_extreme\tCG_pvalue\tCG_zscore\tCG_5perc\tCG_median\tCG_95perc\tCG_mean\tCG_std\n')

	# for each tissue that has a causal gene
	for tissue in CG_tissues:

		print(tissue)
		CSN_nodes = set(CSNs[tissue].vs['name'])
		CGs_in_CSN = tissue_CG_dict[tissue].intersection(CSN_nodes)
		CG_TiP_PPIs = set()
		CGs_with_TiP_PPIs = set()

		TiP_genes_tissue = set(TiPmatrix.loc[TiPmatrix[tissue]>=TiP_cutoff,].index.tolist())
		TiP_genes_CSN = TiP_genes_tissue.intersection(CSN_nodes)

		for CG in CGs_in_CSN:
			neighbors = CSNs[tissue].neighbors(CG)
			n_gene_ids = set([CSNs[tissue].vs[n]['name'] for n in neighbors])
			int_TiPs = n_gene_ids.intersection(TiP_genes_CSN)
			for TiP in int_TiPs:
				if TiP < CG:
					pair = (TiP,CG)
				else:
					pair = (CG,TiP)
				CG_TiP_PPIs.add(pair)
			if len(int_TiPs) > 0:
				CGs_with_TiP_PPIs.add(CG)

		rand_CG_TiP_PPIs_list = []
		rand_CGs_with_TiP_PPIs_list = []

		if len(CGs_in_CSN) > 0:

			# for every random network
			for r in range(num_rand):
				print(r)

				rand_file = config.analysis_path + network_name + '/random_CSN_networks/' + tissue + '/' + network_name + '_rand_network_' + str(r) + '.txt'
				rg = random_graphs.get_real_network(rand_file,None)
				rg_nodes = set(rg.vs['name'])

				rand_CG_TiP_PPIs = set()
				rand_CGs_with_TiP_PPIs = set()

				TiP_genes_rg = TiP_genes_tissue.intersection(rg_nodes)

				for CG in CGs_in_CSN:
					neighbors = rg.neighbors(CG)
					n_gene_ids = set([rg.vs[n]['name'] for n in neighbors])
					int_TiPs = n_gene_ids.intersection(TiP_genes_rg)
					for TiP in int_TiPs:
						if TiP < CG:
							pair = (TiP,CG)
						else:
							pair = (CG,TiP)
						rand_CG_TiP_PPIs.add(pair)
					if len(int_TiPs) > 0:
						rand_CGs_with_TiP_PPIs.add(CG)

				rand_CG_TiP_PPIs_list.append(len(rand_CG_TiP_PPIs))
				rand_CGs_with_TiP_PPIs_list.append(len(rand_CGs_with_TiP_PPIs))

		if len(CGs_in_CSN) == 0:
			target.write(tissue + '\t' + str(len(TiP_genes_CSN)) + '\t' + str(len(CGs_in_CSN)) + '\t' + \
						 str(num_rand) + '\t' + '\t'.join(['NA' for i in range(18)]) + '\n')
		else:
			if numpy.std(rand_CG_TiP_PPIs_list) == 0.0:
				PPI_zscore = 0.0
			else:
				PPI_zscore = (len(CG_TiP_PPIs) - numpy.mean(rand_CG_TiP_PPIs_list))/numpy.std(rand_CG_TiP_PPIs_list)
			if numpy.std(rand_CGs_with_TiP_PPIs_list) == 0.0:
				CG_zscore = 0.0
			else:
				CG_zscore = (len(CGs_with_TiP_PPIs) - numpy.mean(rand_CGs_with_TiP_PPIs_list))/numpy.std(rand_CGs_with_TiP_PPIs_list)
			PPI_num_occ = len(list(filter(lambda v: v >= len(CG_TiP_PPIs),rand_CG_TiP_PPIs_list)))
			PPI_pvalue = PPI_num_occ/float(len(rand_CG_TiP_PPIs_list))
			CG_num_occ = len(list(filter(lambda v: v >= len(CGs_with_TiP_PPIs),rand_CGs_with_TiP_PPIs_list)))
			CG_pvalue = CG_num_occ/float(len(rand_CGs_with_TiP_PPIs_list))

			target.write(tissue + '\t' + str(len(TiP_genes_CSN)) + '\t' + str(len(CGs_in_CSN)) + '\t' + \
						 str(num_rand) + '\t' + str(len(CG_TiP_PPIs)) + '\t' + \
						 str(PPI_num_occ) + '\t' + str(PPI_pvalue) + '\t' + str(PPI_zscore) + '\t' + \
						 str(numpy.percentile(rand_CG_TiP_PPIs_list,5)) + '\t' + \
						 str(numpy.median(rand_CG_TiP_PPIs_list)) + '\t' + \
						 str(numpy.percentile(rand_CG_TiP_PPIs_list,95)) + '\t' + \
						 str(numpy.mean(rand_CG_TiP_PPIs_list)) + '\t' + str(numpy.std(rand_CG_TiP_PPIs_list)) + '\t' + \
						 str(len(CGs_with_TiP_PPIs)) + '\t' + \
						 str(CG_num_occ) + '\t' + str(CG_pvalue) + '\t' + str(CG_zscore) + '\t' + \
						 str(numpy.percentile(rand_CGs_with_TiP_PPIs_list,5)) + '\t' + \
						 str(numpy.median(rand_CGs_with_TiP_PPIs_list)) + '\t' + \
						 str(numpy.percentile(rand_CGs_with_TiP_PPIs_list,95)) + '\t' + \
						 str(numpy.mean(rand_CGs_with_TiP_PPIs_list)) + '\t' + str(numpy.std(rand_CGs_with_TiP_PPIs_list)) + '\n')

	target.close()
