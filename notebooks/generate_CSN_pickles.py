# script to save as pickles various python objects of the networks integrated with GTEx data

import sys
import os
import pickle

import config
import database_utils
import condition_network_classes_v6

if __name__ == '__main__':

	if len(sys.argv) < 3:
		print("""usage:
		1. argument: network name for network to be taken as input or None if all networks are to be run
		Accepted network names:
		2. argument: True for including testis, False for excluding testis""")
		print(config.network_names)
		sys.exit()
	else:
		network_name = sys.argv[1]
		with_testis = sys.argv[2]

	connect = database_utils.get_connection()
	cursor = connect.cursor()

	if network_name == 'None':
		source_networks = config.network_names
	else:
		source_networks = [network_name]


	if with_testis == 'True':
		expr_file = config.expr_file
		suffix = ''
	else:
		expr_file = config.NT_expr_file
		suffix = '_no_testis'

	GTEx = condition_network_classes_v6.Transcriptome(expr_file,config.gene_expr_cutoff)
	GTEx.load_expression_data_ENSG()

	sn_graphs = []
	en_graphs = []
	for network_name in source_networks:

		SN = condition_network_classes_v6.Network(network_name)
		SN.sn = config.load_network(network_name)
		SN.get_expression_annotations_for_nodes(GTEx)
		edge_expr_file = config.analysis_path + network_name + '/GTExv6_' + network_name + '_CSNs.txt'
		SN.get_expression_annotations_for_edges(edge_expr_file)
		SN.get_condition_specific_networks_frac_samples(config.edge_expr_cutoff)
		SN.get_collapsed_CSN_network(GTEx)

		collapsed_CSN_pickle_file = config.analysis_path + network_name + '/GTExv6_collapsed_CSN_' + network_name + suffix + '.pickle'
		target = open(collapsed_CSN_pickle_file,'wb')
		pickle.dump(SN.collapsed_CSN,target)
		target.close()
