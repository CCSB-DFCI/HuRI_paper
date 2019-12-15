# script that generates tissue-filtered networks, DB tables to store them and
# loading functions

import sys
import igraph, pandas
import ccsblib.huri
import database_utils
from custom_settings import *

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
def load_network(network_name,cache='read',no_homodimer_ppis=True):

	if network_name == 'HI-I-05':
		snw = ccsblib.huri.load_nw_hi_union(id_type='ensembl_gene_id',fmt='igraph',cache=cache)
		hi_i_05_edges = snw.es.select(lambda e: e['in_Rual'])
		hi_i_05_edges_ensg_id = [(snw.vs[e.source]['name'],snw.vs[e.target]['name']) for e in hi_i_05_edges]
		nw = igraph.Graph.TupleList(hi_i_05_edges_ensg_id)
		for attribute in hi_i_05_edges.attributes():
			nw.es[attribute] = hi_i_05_edges[attribute]
	elif network_name == 'HI-II-14':
		snw = ccsblib.huri.load_nw_hi_union(id_type='ensembl_gene_id',fmt='igraph',cache=cache)
		hi_ii_14_edges = snw.es.select(lambda e: e['in_HI_II_14_screen_1'] or e['in_HI_II_14_screen_2'])
		hi_ii_14_edges_ensg_id = [(snw.vs[e.source]['name'],snw.vs[e.target]['name']) for e in hi_ii_14_edges]
		nw = igraph.Graph.TupleList(hi_ii_14_edges_ensg_id)
		for attribute in hi_ii_14_edges.attributes():
			nw.es[attribute] = hi_ii_14_edges[attribute]
	elif network_name == 'HI-III':
		nw = ccsblib.huri.load_nw_hi_iii(id_type='ensembl_gene_id',fmt='igraph',cache=cache)
	elif network_name == 'HI-union':
		nw = ccsblib.huri.load_nw_hi_union_for_paper(id_type='ensembl_gene_id',fmt='igraph',cache=cache)
	elif network_name == 'Lit-BM-17':
		nw = ccsblib.huri.load_nw_lit_bm_17(id_type='ensembl_gene_id',fmt='igraph',cache=cache)
	elif network_name == 'BioPlex':
		nw = ccsblib.huri.load_nw_bioplex(id_type='ensembl_gene_id',fmt='igraph',cache=cache)
	elif network_name == 'CoFrac':
		nw = ccsblib.huri.load_nw_cofrac(id_type='ensembl_gene_id',fmt='igraph',cache=cache)
	elif network_name == 'QUBIC':
		nw = ccsblib.huri.load_nw_qubic(id_type='ensembl_gene_id',fmt='igraph',cache=cache)
	else:
		print('Wrong network name.')
		sys.exit()

	if no_homodimer_ppis:
		nw = remove_homodimer_ppis(nw)

	return nw


# function to generate TSN
def generate_TSN(network_name,cursor):

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
	query = """select tissue_names from {}.GTEx_tissue_names order by tissue_names asc""".format(DB_NAME)
	cursor.execute(query)
	tissues = [row[0] for row in cursor]

	# for every tissue:
	for tissue in tissues:
		# load expression data
		df = pandas.read_table('../external_data/GTEx_sample_data/' + tissue + '.txt',index_col=0)
		df.rename(mapper=(lambda x: x.split('.')[0]),axis='index',inplace=True)
		df = df.applymap(lambda x: 1 if x > GTEX_EXPR_CUTOFF else 0)
		# calculate fraction samples where edge is co-expressed for all edges
		# add to network dataframe
		nw_df[tissue] = nw_df.apply(lambda x: get_frac_samples(x,df),axis=1)

	# write out to DB table
	DB_TABLE = network_name.replace('-','_') + '_TSNs'
	query = """delete from {}.{}""".format(DB_NAME,DB_TABLE)
	cursor.execute(query)
	for i in range(nw_df.shape[0]):
	    query = 'insert into {}.{} values ('.format(DB_NAME,DB_TABLE) + ','.join(['%s' for c in range(nw_df.shape[1])]) + ')'
	    cursor.execute(query,tuple(nw_df.iloc[i,].tolist()))


if __name__ == '__main__':

	connect = database_utils.get_connection()
	cursor = connect.cursor()

	for network_name in NETWORK_NAMES:
		print(network_name)
		generate_TSN(network_name,cursor)
