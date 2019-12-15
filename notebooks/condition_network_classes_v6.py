
import igraph
import MySQLdb
import random
import numpy

class Network:

	def __init__(self,name):

		self.sn = igraph.Graph(n=0,vertex_attrs={'name':[]})
		self.sn['name'] = name
		self.conditions = set()
		self.CSNs = {}
		self.sb_CSNs = {}
		self.edge_condition_dict = {}
		self.collapsed_CSN = igraph.Graph(n=0,vertex_attrs={'name':[]})
		self.space = None
		self.space_name = None


	# function that adds edges and the respective nodes from the given
	# source networks
	def get_source_network(self,source_networks,cursor):

		# stores all geneIDs
		nodes = set()
		# stores for every edge a set of source names
		edge_dict = {}

		for tup in source_networks:
			source_name = tup[0]
			source = tup[1]
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
					if edge not in edge_dict:
						edge_dict[edge] = set()
					edge_dict[edge].add(source_name)
			# the source is a file
			else:
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
					if edge not in edge_dict:
						edge_dict[edge] = set()
					edge_dict[edge].add(source_name)
		for edge,origins in edge_dict.items():
			nodes.add(edge[0])
			nodes.add(edge[1])
		# add all nodes to the network
		self.sn.add_vertices(list(nodes))
		# add all edges to the network together with the origins attribute
		items = edge_dict.items()
		edges = [tup[0] for tup in items]
		ori_lists = [list(tup[1]) for tup in items]
		self.sn.add_edges(edges)
		self.sn.es['origins'] = ori_lists


	# function that stores expression information for all nodes from given input file
	def get_expression_annotations_for_nodes(self,transcriptome):

		node_names = self.sn.vs['name']
		expr_dict = {}

		for condition in transcriptome.condition_expr_dict.keys():
			expr_dict[condition] = [0.0 for i in range(0,len(node_names))]

		for geneID,condition_dict in transcriptome.gene_expr_dict.items():
			if geneID in node_names:
				index = node_names.index(geneID)
				for condition, value in condition_dict.items():
					expr_dict[condition][index] = float(value)

		for condition,values in expr_dict.items():
			self.sn.vs[condition] = values
			self.conditions.add(condition)


	# function that stores in every edge a list with all the fractions of samples per condition in
	# which the edge was observed
	def get_expression_annotations_for_edges(self,infile):

		file1 = open(infile,'r')
		entries = file1.readlines()
		file1.close()

		num_edges = len(self.sn.es)

		expr_dict = {}

		if len(self.conditions) == 0:
			print('Warning, no conditions known.')

		for condition in self.conditions:
			expr_dict[condition] = [0.0 for i in range(0,num_edges)]

		input_conditions = entries[0][:-1].split('\t')[2:]

		for line in entries[1:]:
			tab_list = str.split(line[:-1],'\t')
			geneA = tab_list[0]
			indexA = self.sn.vs.find(geneA).index
			geneB = tab_list[1]
			indexB = self.sn.vs.find(geneB).index
			edge_index = self.sn.get_eid(indexA,indexB)
			for j,value in enumerate(tab_list[2:]):
				if input_conditions[j] in self.conditions:
					expr_dict[input_conditions[j]][edge_index] = float(value)

		for condition,values in expr_dict.items():
			self.sn.es[condition] = values


	# function that derives condition-specific networks based on in how many samples per condition
	# the co-expressed genes were observed
	def get_condition_specific_networks_frac_samples(self,cutoff):

		for condition in self.conditions:
			csn_edges = self.sn.es.select(lambda e: e[condition] > cutoff)
			self.CSNs[condition] = self.sn.subgraph_edges(csn_edges)


	# function that creates a dictionary listing for every edge between co-expressed genes from the
	# source network in which tissues this edge has been seen
	def get_edge_condition_dict(self,Transcriptome):

		for condition,csn in self.CSNs.items():
			for edge_int in csn.get_edgelist():
				edge_id_a = csn.vs[edge_int[0]]['name']
				edge_id_b = csn.vs[edge_int[1]]['name']
				if edge_id_a < edge_id_b:
					edge_id = (edge_id_a,edge_id_b)
				else:
					edge_id = (edge_id_b,edge_id_a)
				if edge_id not in self.edge_condition_dict:
					self.edge_condition_dict[edge_id] = set()
				self.edge_condition_dict[edge_id].add(condition)


	# function that builds a subgraph consisting of all edges that are between
	# co-expressed genes
	def get_collapsed_CSN_network(self,Transcriptome=None):

		if len(self.edge_condition_dict) == 0:
			self.get_edge_condition_dict(Transcriptome)

		# get all nodes
		geneIDs = set()
		edges = self.edge_condition_dict.keys()
		for tup in edges:
			geneIDs.add(tup[0])
			geneIDs.add(tup[1])
		# add all nodes to the network
		self.collapsed_CSN.add_vertices(list(geneIDs))
		# add all edges to the network
		self.collapsed_CSN.add_edges(edges)


class Transcriptome:

	def __init__(self,infile,cutoff):

		self.expr_file = infile
		self.cutoff = cutoff
		self.gene_expr_dict = {}
		self.condition_expr_dict = {}
		self.gene_sets = {}
		self.gene_sets_per_condition = {}
		self.TiP_dict = {}
		self.TiP_gene_expr_dict = {}


	def load_expression_data_ENSG(self):

		file1 = open(self.expr_file,'r')
		entries = file1.readlines()
		file1.close()
		conditions = entries[0][:-1].split('\t')[1:]
		for condition in conditions:
			self.condition_expr_dict[condition] = {}
		for line in entries[1:]:
			tab_list = str.split(line[:-1],'\t')
			geneID = tab_list[0]
			values = [float(i) for i in tab_list[1:]]
			if max(values) > self.cutoff:
				self.gene_expr_dict[geneID] = {}
				for i,value in enumerate(values):
					if value <= self.cutoff:
						self.gene_expr_dict[geneID][conditions[i]] = 0.0
						self.condition_expr_dict[conditions[i]][geneID] = 0.0
					else:
						self.gene_expr_dict[geneID][conditions[i]] = float(value)
						self.condition_expr_dict[conditions[i]][geneID] = float(value)
