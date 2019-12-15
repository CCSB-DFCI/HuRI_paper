# script to create files for drawing networks around TiP genes

import os, sys, pickle
import pandas, igraph, numpy
import database_utils, config
import TiP_network_closeness

def get_gene_symbol_descr_info(cursor):

	query = """select a.ensembl_gene_id_short,a.symbol,b.DESCRIPTION from
	        horfeome_annotation_gencode27.gencode2entrez a
	        left join
	        entrez_gene.gene_info_human_20170914 b
	        on a.entrez_gene_id=b.GENE_ID"""
	cursor.execute(query)
	gene_name_dict = {}
	for row in cursor:
	    gene_id = row[0]
	    symbol = row[1]
	    descr = row[2]
	    gene_name_dict[gene_id] = (symbol,descr)

	return gene_name_dict


# function that writes out the nodes and their attributes from the given graph object in cytoscape readable format
def write_node_attributes(graph,outfile):

	target = open(outfile,'w')
	node_attributes = graph.vertex_attributes()
	node_attributes.sort()
	target.write('geneID')
	for attr in node_attributes:
		if attr != 'name':
			target.write('\t' + attr)
	target.write('\n')

	for node in graph.vs:
		target.write(node['name'])
		for attr in node_attributes:
			if attr != 'name':
				if isinstance(node[attr],list):
					outstring = ''
					for value in node[attr]:
						outstring = outstring + str(value) + '|'
					outstring = outstring[:-1]
				else:
					outstring = str(node[attr])
				target.write('\t' + outstring)
		target.write('\n')

	target.close()


# function that writes out the edges and their attributes from the given graph object in cytoscape readable format
def write_edge_attributes(graph,outfile):

	target = open(outfile,'w')
	edge_attributes = graph.edge_attributes()
	edge_attributes.sort()
	target.write('geneID_A\tgeneID_B')
	for attr in edge_attributes:
		target.write('\t' + attr)
	target.write('\n')

	for edge in graph.es:
		source_node = graph.vs[edge.source]['name']
		target_node = graph.vs[edge.target]['name']
		if source_node < target_node:
			target.write(str(source_node) + '\t' + str(target_node))
		else:
			target.write(str(target_node) + '\t' + str(source_node))
		for attr in edge_attributes:
			if isinstance(edge[attr],list):
				outstring = ''
				for value in edge[attr]:
					outstring = outstring + str(value) + '|'
				outstring = outstring[:-1]
			else:
				outstring = str(edge[attr])
			target.write('\t' + outstring)
		target.write('\n')

	target.close()


if __name__ == '__main__':

	network_name = sys.argv[1]
	tissues = sys.argv[2].split(',')
	output_folder = sys.argv[3]

	connect = database_utils.get_connection()
	cursor = connect.cursor()

	file1 = open(config.NT_cond_file,'r')
	all_tissues = file1.readlines()
	file1.close()
	all_tissues = [c[:-1] for c in all_tissues]

	if len(tissues) == 1 and tissues[0] == 'all':
		tissues = all_tissues

	print('load gene expression data')
	GTEx = pandas.read_table(config.NT_expr_file,index_col=0)
	TiPmatrix = pandas.read_table(config.NT_pref_expr_file,index_col=0)
	GTEx_tissues = TiPmatrix.columns
	for tissue in GTEx_tissues:
		TiPmatrix.loc[TiPmatrix.index.isin(GTEx.index[GTEx[tissue]<=config.gene_expr_cutoff]),[tissue]] = numpy.NaN

	print('Read in gene annotation data')
	gene_symbol_descr_dict = get_gene_symbol_descr_info(cursor)

	tissue_graphs = {}

	for tissue in tissues:

		print(tissue)

		tisTiPgraph = igraph.Graph()
		tisTiPgraph.vs['name'] = []

		CSN = TiP_network_closeness.get_CSN_for_tissue_and_source_network(network_name,tissue)
		CSN_nodes = set(CSN.vs['name'])

		# get the TiP genes for this tissue
		tissue_TiPs = set(TiPmatrix.loc[TiPmatrix[tissue]>=2,].index.tolist())

		# build a graph linking all TiP genes with their neighbors
		# get the neighbors for the TiP genes and add the corresponding nodes and edges to the network
		TiP_genes_CSN = tissue_TiPs.intersection(CSN_nodes)
		tisGraph_nodes = set()
		tisGraph_nodes = tisGraph_nodes.union(TiP_genes_CSN)
		tisGraph_edges = set()
		for TiPgene in TiP_genes_CSN:
			n_list = CSN.neighbors(TiPgene)
			n_gene_ids = [CSN.vs[n]['name'] for n in n_list]
			tisGraph_nodes = tisGraph_nodes.union(set(n_gene_ids))
			n_edges = [tuple(sorted([TiPgene,n_gene_id])) for n_gene_id in n_gene_ids]
			tisGraph_edges = tisGraph_edges.union(set(n_edges))
		tisTiPgraph.add_vertices(list(tisGraph_nodes))
		tisTiPgraph.add_edges(list(tisGraph_edges))

		# add TiP values
		tisTiPgraph.vs['TiP_value'] = [0.0 for n in range(len(tisTiPgraph.vs))]
		for n in tisTiPgraph.vs:
			n['TiP_value'] = TiPmatrix.loc[TiPmatrix.index==n['name'],[tissue]].values[0][0]
		discr_TiP_values = [(0 if n['TiP_value'] < 2 else n['TiP_value']) for n in tisTiPgraph.vs]
		discr_TiP_values = [(2 if (value >= 2 and value < 3) else value) for value in discr_TiP_values]
		discr_TiP_values = [(3 if (value >= 3 and value < 9) else value) for value in discr_TiP_values]
		discr_TiP_values = [(9 if value >= 9 else value) for value in discr_TiP_values]
		tisTiPgraph.vs['TiP_value_discrete'] = discr_TiP_values

		# add info which edge links two TiP genes
		values = []
		for edge in tisTiPgraph.es:
			TiPvalue_a = tisTiPgraph.vs[edge.source]['TiP_value']
			TiPvalue_b = tisTiPgraph.vs[edge.target]['TiP_value']
			if TiPvalue_a >= 2 and TiPvalue_b >= 2:
				values.append(1)
			else:
				values.append(0)
		tisTiPgraph.es['TiP_TiP_link'] = values

		# gene symbol and description information
		tisTiPgraph.vs['symbol'] = [None for n in range(len(tisTiPgraph.vs))]
		tisTiPgraph.vs['descr'] = [None for n in range(len(tisTiPgraph.vs))]
		for n in range(len(tisTiPgraph.vs)):
			geneID = tisTiPgraph.vs[n]['name']
			if geneID in gene_symbol_descr_dict:
				tisTiPgraph.vs[n]['symbol'] = gene_symbol_descr_dict[geneID][0]
				tisTiPgraph.vs[n]['descr'] = gene_symbol_descr_dict[geneID][1]

		# add gene expression values
		tisTiPgraph.vs['expr'] = [None for n in range(len(tisTiPgraph.vs))]
		for n in range(len(tisTiPgraph.vs)):
			geneID = tisTiPgraph.vs[n]['name']
			tisTiPgraph.vs[n]['expr'] = GTEx.loc[GTEx.index==geneID,[tissue]].values[0][0]

		node_outfile = config.analysis_path + network_name + '/' + output_folder + '/' + network_name + '_' + tissue + '.node_attributes.txt'
		write_node_attributes(tisTiPgraph,node_outfile)
		edge_outfile = config.analysis_path + network_name + '/' + output_folder + '/' + network_name + '_' + tissue + '.edge_attributes.txt'
		write_edge_attributes(tisTiPgraph,edge_outfile)

		tissue_graphs[tissue] = tisTiPgraph

		print(tissue, 'nodes:', len(tisTiPgraph.vs), 'edges:', len(tisTiPgraph.es), '#tissue TiPs:', len(tissue_TiPs))

	target = open(config.analysis_path + network_name + '/' + output_folder + '/' + output_folder + '.pickle','bw')
	pickle.dump(tissue_graphs,target)
	target.close()
