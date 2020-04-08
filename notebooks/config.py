# config file

import os, sys
import igraph

from utils import (load_nw_hi_iii,
					load_nw_hi_union,
					load_nw_lit_bm_17,
					load_nw_qubic,
					load_nw_bioplex,
					load_nw_cofrac)


rand_home_dir = '../data/katjas_data'
path = '../data/katjas_data/GTEx/'
data_path = path + 'expression_datasets/'
analysis_path = path + 'analysis_GTExV6_2018_frozen/'
global_path = analysis_path + 'global/'
rand_dir = '../data/katjas_data/GTEx/analysis_GTExV6_2018_frozen/'
cond_file = data_path + 'GTExv6_JPnorm_subtypes_no_cell.txt'
NT_cond_file = data_path + 'GTExv6_JPnorm_subtypes_no_cell_no_testis.txt'
expr_file = data_path + 'GTExv6_JPnorm_gene_expression_per_subtype_ENSG-PC.txt'
NT_expr_file = data_path + 'GTExv6_JPnorm_gene_expression_per_subtype_ENSG-PC_no_testis.txt'
pref_expr_file = data_path + 'GTExv6_tissue_specificity_data_no_cutoff_PC_no_cell.txt'
NT_pref_expr_file = data_path + 'GTExv6_tissue_specificity_data_no_cutoff_PC_no_cell_no_testis.txt'
gene_expr_cutoff = 5
edge_expr_cutoff = 0.5
# network_names = ['HI-I-05','HI-II-14','HI-III-16','HI-III-16-2ev','HI-union-16','HI-III','HI-III-2ev','HI-union',\
# 				 'Lit-BM-17','BioPlex','CoFrac','QUBIC']
network_names = ['HI-III','HI-union','Lit-BM-17','BioPlex','QUBIC','CoFrac']


def remove_homodimer_ppis(g):

	edges_to_delete = []
	for edge in g.es:
		if edge.source == edge.target:
			edges_to_delete.append(edge)
	g.delete_edges(edges_to_delete)

	return g


def ev2_edge(edge):
    attr = edge.attributes().items()
    screens = sum([t[1] for t in filter(lambda t: t[0].find('in_screen') > -1,attr)])
    assays = sum([t[1] for t in filter(lambda t: t[0].find('in_assay') > -1,attr)])
    if screens > 1 or assays > 1:
        return True
    else:
        return False


def get_2ev_network(g):

	ev2_edges = g.es.select(lambda e: ev2_edge(e))
	ev2_edges_ensg_id = [(g.vs[e.source]['name'],g.vs[e.target]['name']) for e in ev2_edges]
	g_ev2 = igraph.Graph.TupleList(ev2_edges_ensg_id)
	for attribute in ev2_edges.attributes():
		g_ev2.es[attribute] = ev2_edges[attribute]

	return g_ev2


def load_networks(network_names=network_names,cache='read',no_homodimer_ppis=True):

	networks = {}
	for network_name in network_names:
		networks[network_name] = load_network(network_name,cache=cache,no_homodimer_ppis=no_homodimer_ppis)

	return networks


def load_network(network_name,cache='read',no_homodimer_ppis=True):

	if network_name == 'HI-I-05':
		snw = load_nw_hi_union(id_type='ensembl_gene_id',fmt='igraph')
		hi_i_05_edges = snw.es.select(lambda e: e['in_Rual'])
		hi_i_05_edges_ensg_id = [(snw.vs[e.source]['name'],snw.vs[e.target]['name']) for e in hi_i_05_edges]
		nw = igraph.Graph.TupleList(hi_i_05_edges_ensg_id)
		for attribute in hi_i_05_edges.attributes():
			nw.es[attribute] = hi_i_05_edges[attribute]
	elif network_name == 'HI-II-14':
		snw = load_nw_hi_union(id_type='ensembl_gene_id',fmt='igraph')
		hi_ii_14_edges = snw.es.select(lambda e: e['in_HI_II_14_screen_1'] or e['in_HI_II_14_screen_2'])
		hi_ii_14_edges_ensg_id = [(snw.vs[e.source]['name'],snw.vs[e.target]['name']) for e in hi_ii_14_edges]
		nw = igraph.Graph.TupleList(hi_ii_14_edges_ensg_id)
		for attribute in hi_ii_14_edges.attributes():
			nw.es[attribute] = hi_ii_14_edges[attribute]
	elif network_name == 'HI-III':
		nw = load_nw_hi_iii(id_type='ensembl_gene_id',fmt='igraph')
	elif network_name == 'HI-III-2ev':
		nw = load_nw_hi_iii(id_type='ensembl_gene_id',fmt='igraph')
		nw = get_2ev_network(nw)
	elif network_name == 'HI-union':
		nw = load_nw_hi_union(id_type='ensembl_gene_id',fmt='igraph')
	elif network_name == 'Lit-BM-17':
		nw = load_nw_lit_bm_17(id_type='ensembl_gene_id',fmt='igraph')
	elif network_name == 'BioPlex':
		nw = load_nw_bioplex(id_type='ensembl_gene_id',fmt='igraph')
	elif network_name == 'CoFrac':
		nw = load_nw_cofrac(id_type='ensembl_gene_id',fmt='igraph')
	elif network_name == 'QUBIC':
		nw = load_nw_qubic(id_type='ensembl_gene_id',fmt='igraph')
	elif network_name == 'HI_BP':
		nw = load_nw_hi_iii(id_type='ensembl_gene_id',fmt='igraph')
		bp = load_nw_bioplex(id_type='ensembl_gene_id',fmt='igraph')
		new_nodes = set(bp.vs['name']).difference(set(nw.vs['name']))
		nw.add_vertices(list(new_nodes))
		new_edges = set()
		for edge in bp.es:
			gene_a = bp.vs[edge.source]['name']
			gene_b = bp.vs[edge.target]['name']
			if not nw.are_connected(gene_a,gene_b):
				new_edges.add((gene_a,gene_b))
		nw.add_edges(list(new_edges))
	else:
		print('Wrong network name.')
		sys.exit()

	if no_homodimer_ppis:
		nw = remove_homodimer_ppis(nw)

	return nw


def write_networks(networks):

	subfolders = os.listdir(analysis_path)
	for network_name,g in networks.items():
		if network_name not in subfolders:
			os.mkdir(analysis_path + network_name)
		outfile = analysis_path + network_name + '/' + network_name + '.tsv'
		write_network(g,outfile)


def write_network(g,outfile):

		target = open(outfile,'w')
		target.write('ensembl_gene_id_a\tensembl_gene_id_b\n')
		for edge in g.es:
			gene_a = g.vs[edge.source]['name']
			gene_b = g.vs[edge.target]['name']
			if gene_a < gene_b:
				target.write(gene_a + '\t' + gene_b + '\n')
			else:
				target.write(gene_b + '\t' + gene_a + '\n')
		target.close()


def get_CSNs_for_source_network(name_source_network,analysis_path=analysis_path,edge_expr_cutoff=edge_expr_cutoff):

	file1 = open(analysis_path + name_source_network + '/GTExv6_' + name_source_network + '_CSNs.txt','r')
	entries = file1.readlines()
	file1.close()
	CSNs = {}
	tissues = entries[0][:-1].split('\t')[2:]
	for t,tissue in enumerate(tissues):
		CSNs[tissue] = get_CSNs_for_tissue_and_source_network(entries,t)

	return CSNs


def get_CSNs_for_tissue_and_source_network(CSN_file_entries,tis_index,edge_expr_cutoff=edge_expr_cutoff):

	edge_list = []
	for line in CSN_file_entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		gene_a = tab_list[0]
		gene_b = tab_list[1]
		frac_samples = float(tab_list[tis_index+2])
		if frac_samples > edge_expr_cutoff:
			edge_list.append((gene_a,gene_b))
	nodes = set()
	for edge in edge_list:
		nodes.add(edge[0])
		nodes.add(edge[1])
	CSN = igraph.Graph()
	CSN.add_vertices(list(nodes))
	CSN.add_edges(edge_list)

	return CSN
