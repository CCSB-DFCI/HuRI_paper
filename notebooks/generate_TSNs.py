# script that generates tissue-filtered networks, DB tables to store them and
# loading functions
import os
import sys
import igraph
import pandas
from custom_settings import *

from utils import (load_nw_hi_iii,
					load_nw_hi_union,
					load_nw_lit_bm_17,
					load_nw_qubic,
					load_nw_bioplex,
					load_nw_cofrac)

# gets igraph object and deletes all homodimer edges from it
def remove_homodimer_ppis(g):

	edges_to_delete = []
	for edge in g.es:
		if edge.source == edge.target:
			edges_to_delete.append(edge)
	g.delete_edges(edges_to_delete)

	orphan_nodes = g.vs.select(_degree_eq=0)
	g.delete_vertices(orphan_nodes)

	return g

# function to load source networks with homodimers removed
def load_network(network_name, no_homodimer_ppis=True):

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
	elif network_name == 'HI-union':
		nw = load_nw_hi_union_for_paper(id_type='ensembl_gene_id',fmt='igraph')
	elif network_name == 'Lit-BM-17':
		nw = load_nw_lit_bm_17(id_type='ensembl_gene_id',fmt='igraph')
	elif network_name == 'BioPlex':
		nw = load_nw_bioplex(id_type='ensembl_gene_id',fmt='igraph')
	elif network_name == 'CoFrac':
		nw = load_nw_cofrac(id_type='ensembl_gene_id',fmt='igraph')
	elif network_name == 'QUBIC':
		nw = load_nw_qubic(id_type='ensembl_gene_id',fmt='igraph')
	else:
		print('Wrong network name.')
		sys.exit()

	if no_homodimer_ppis:
		nw = remove_homodimer_ppis(nw)

	return nw


# function to generate TSN
def generate_TSN(network_name):

	# load source network into dataframe
	nw = load_network(network_name)
	geneAs = []
	geneBs = []
	for edge in nw.es:
		geneA = nw.vs[edge.source]['name']
		geneB = nw.vs[edge.target]['name']
		if geneA < geneB:
			geneAs.append(geneA)
			geneBs.append(geneB)
		else:
			geneAs.append(geneB)
			geneBs.append(geneA)
	nw_df = pandas.DataFrame({'gene_a':geneAs,'gene_b':geneBs})

	def get_frac_samples(x,df):

		if x['gene_a'] in df.index and x['gene_b'] in df.index:
			gene_a_values = df.loc[x['gene_a'],].tolist()
			gene_b_values = df.loc[x['gene_b'],].tolist()
			values = list(map(sum,zip(gene_a_values,gene_b_values)))
			num_expr = len(list(filter(lambda v: v == 2,values)))
			frac = num_expr/float(len(values))
		else:
			frac = 0
		return frac

	# load tissue names
	tissues = pandas.read_csv('../data/processed/Supplementary Table 27.txt', sep='\t').set_index('Ensembl_gene_id').columns

	# for every tissue:
	for tissue in tissues:
		# load expression data
		df = pandas.read_table('../external_data/GTEx_sample_data/' + tissue + '.txt',index_col=0)
		df.rename(mapper=(lambda x: x.split('.')[0]),axis='index',inplace=True)
		df = df.applymap(lambda x: 1 if x > GTEX_EXPR_CUTOFF else 0)
		# calculate fraction samples where edge is co-expressed for all edges
		# add to network dataframe
		nw_df[tissue] = nw_df.apply(lambda x: get_frac_samples(x,df),axis=1)

	# write out to file
	TABLE = network_name.replace('-','_') + '_TSNs'
	out_dir = '../data/processed/tissue_specific_networks'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	nw_df.to_csv(os.path.join(out_dir, TABLE + '.tsv'), sep='\t')

if __name__ == '__main__':
	for network_name in NETWORK_NAMES:
		print(network_name)
		generate_TSN(network_name)
