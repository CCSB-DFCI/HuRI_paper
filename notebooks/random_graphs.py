
# script to generate random graphs

import sys
import os
import pickle
import MySQLdb
import igraph
import database_utils
import config
home_dir = os.environ['HOME']

def get_real_network(source,cursor):

	# stores all geneIDs
	nodes = set()
	edge_list = []

	# if the source is a binary pickle file that contains the igraph object of the graph, read in that object
	if source[-6:] == 'pickle':
		file1 = open(source,'rb')
		g = pickle.load(file1)
	else:
		# if the source is a database query
		if source[:6] == 'select':
			cursor.execute(source)
			rows = cursor.fetchall()
			for row in rows:
				gene_a = str(row[0])
				gene_b = str(row[1])
				if gene_a < gene_b:
					edge = (gene_a,gene_b)
				else:
					edge = (gene_b,gene_a)
				edge_list.append(edge)
		# the source is a simple text file with edges
		elif source[-3:] == 'txt' or source[-3:] == 'tsv':
			file1 = open(source,'r')
			entries = file1.readlines()
			file1.close()
			for line in entries[1:]:
				tab_list = str.split(line[:-1],'\t')
				gene_a = tab_list[0]
				gene_b = tab_list[1]
				if gene_a < gene_b:
					edge = (gene_a,gene_b)
				else:
					edge = (gene_b,gene_a)
				edge_list.append(edge)
		else:
			print('unknown source type')
			sys.exit()

		for edge in edge_list:
			nodes.add(edge[0])
			nodes.add(edge[1])

		g = igraph.Graph()
		# add all nodes to the network
		g.add_vertices(list(nodes))
		# add all edges to the network
		g.add_edges(edge_list)

	return g


def get_random_graph(real_network,outfile):

	degree_seq = real_network.degree()
	rand_graph = igraph.Graph.Degree_Sequence(degree_seq,method='vl')
	rand_graph.vs['name'] = real_network.vs['name']
	rand_edges = rand_graph.get_edgelist()
	rand_edges_mapped = []
	for edge in rand_edges:
		gene_a = real_network.vs[edge[0]]['name']
		gene_b = real_network.vs[edge[1]]['name']
		if gene_a < gene_b:
			rand_edges_mapped.append((gene_a,gene_b))
		else:
			rand_edges_mapped.append((gene_b,gene_a))
	target = open(outfile,'w')
	target.write('gene_a\tgene_b\n')
	for edge in rand_edges_mapped:
		target.write(str(edge[0]) + '\t' + str(edge[1]) + '\n')
	target.close()


def get_random_graphs(source,outpath,source_name,num_graphs,cursor):

	real_network = get_real_network(source,cursor)

	for i in range(0,num_graphs):
		if i%100 == 0:
			print(i)
		outfile = outpath + source_name + '_rand_network_' + str(i) + '.txt'
		get_random_graph(real_network,outfile)


if '__main__' == __name__:

	connect = database_utils.get_connection()
	cursor = connect.cursor()

	network_name = sys.argv[1]
	num_graphs = int(sys.argv[2])
	rand_path = sys.argv[3]
	rand_folder = sys.argv[4]
	tissue_spec = sys.argv[5]

	if not os.path.exists(rand_path + network_name):
		os.mkdir(rand_path + network_name)
	if not os.path.exists(rand_path + network_name + '/' + rand_folder):
		os.mkdir(rand_path + network_name + '/' + rand_folder)

	if tissue_spec == '0':
		# code to generate random graphs for source networks
		network_file = config.analysis_path + network_name + '/' + network_name + '.tsv'
		outpath = rand_path + network_name + '/' + rand_folder + '/'
		get_random_graphs(network_file,outpath,network_name,num_graphs,cursor)

	elif tissue_spec == '1':
		# code to generate random networks for every specific CSN of every original network
		# read in the CSNs as dict
		# go through that dict and generate a 1000 random networks for every tissue
		# create directory that will contain all random networks
		network_path = rand_path + network_name + '/'
		CSNs = config.get_CSNs_for_source_network(network_name)
		for condition,CSN in CSNs.items():
			r_subdirs = os.listdir(network_path + rand_folder)
			if condition not in r_subdirs:
				os.mkdir(network_path + rand_folder + '/' + condition)
			cond_path = network_path + rand_folder + '/' + condition + '/'
			for r in range(0,num_graphs):
				outfile = cond_path + network_name + '_rand_network_' + str(r) + '.txt'
				get_random_graph(CSN,outfile)

	elif tissue_spec == '2':
		# code to generate random networks from the collapsed CSNs
		# starting from a pickle file
		with_testis = sys.argv[6]
		if with_testis == 'True':
			suffix = ''
		else:
			suffix = '_no_testis'

		source = config.analysis_path + network_name + '/GTExv6_collapsed_CSN_' + network_name + suffix + '.pickle'
		outpath = rand_path + network_name + '/' + rand_folder + '/'
		get_random_graphs(source,outpath,network_name,num_graphs,cursor)


	else:
		print('Unknown input option. Exit')
