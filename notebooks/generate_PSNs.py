# script to generate PSNs from real and random networks

import igraph
import sys, os
home_dir = os.environ['HOME']

import random_graphs


def build_PSN(g):

	psn = igraph.Graph(n=0,vertex_attrs={'name':[]},edge_attrs={'RAI':[]})
	psn.add_vertices(g.vs['name'])
	RAIs = []
	JSs = []
	num_shared_partners = []
	edges = []
	d2_nodes = g.vs.select(_degree_ge=2)
	for n1, node1 in enumerate(d2_nodes):
		n1_neighbors = set(g.neighbors(node1))
		for node2 in d2_nodes[n1+1:]:
			n2_neighbors = set(g.neighbors(node2))
			common_neighbors = n1_neighbors.intersection(n2_neighbors)
			if len(common_neighbors) > 0:
				RAIs.append(sum([1.0/g.vs[cn].degree() for cn in common_neighbors]))
				JSs.append(len(common_neighbors)/float(len(n1_neighbors.union(n2_neighbors))))
				num_shared_partners.append(len(common_neighbors))
				edges.append(tuple(sorted([node1,node2])))
	psn.add_edges(edges)
	psn.es['RAI'] = RAIs
	psn.es['JS'] = JSs
	psn.es['num_shared_partners'] = num_shared_partners
	orphan_nodes = psn.vs.select(_degree = 0)
	orphan_nodes_indices = [v.index for v in orphan_nodes]
	psn.delete_vertices(orphan_nodes_indices)

	return psn


if __name__ == '__main__':

	# save every edge between nodes that share at least 1 interaction partner
	# save RAI, Jaccard similarity, # shared interaction partners, number interaction partners a and b
	# save the same for every random HI-III network

	if len(sys.argv) < 7:
		print(\
		"""insufficient number of arguments:
		network name
		generate PSN from real network?
		number of random network start
		number of random network end
		network path
		exclude degree 1 nodes: 1 - yes, 0 - no
		""")
		sys.exit()
	else:
		network_name = sys.argv[1]
		real_PSN = sys.argv[2]
		num_rand_start = int(sys.argv[3])
		num_rand_end = int(sys.argv[4])
		network_path = sys.argv[5]
		deg1_nodes = int(sys.argv[6])

	if real_PSN == 'True':
		network = random_graphs.get_real_network(network_path + network_name + '.txt',None)
		print(len(network.es),len(network.vs))
		if deg1_nodes == 1:
			deg1_vertices = network.vs.select(_degree=1)
			node_indices = [v.index for v in deg1_vertices]
			network.delete_vertices(node_indices)
		print(len(network.es),len(network.vs))

		# get real PSN and write out
		real_PSN = build_PSN(network)
		print(len(real_PSN.es),len(real_PSN.vs))
		target = open(network_path + network_name + '_PSN.txt','w')
		target.write('gene_a\tgene_b\tRAI\tJaccard_similarity\tnum_shared_partners\tnum_partners_a\tnum_partners_b\n')
		for edge in real_PSN.es:
			gene_ids = sorted([real_PSN.vs[edge.source]['name'],real_PSN.vs[edge.target]['name']])
			target.write(gene_ids[0] + '\t' + gene_ids[1] + '\t' + str(edge['RAI']) + '\t' + str(edge['JS']) + '\t' + \
						 str(edge['num_shared_partners']) + '\t' + str(network.degree(gene_ids[0])) + \
						 '\t' + str(network.degree(gene_ids[1])) + '\n')
		target.close()

	rand_path = network_path + network_name + '_PSN/'
	if not os.path.exists(rand_path):
		os.mkdir(rand_path)

	for r in range(num_rand_start,num_rand_end+1):

		print(r)

		rand_file = network_path + network_name + '/' + network_name + '_rand_network_' + str(r) + '.txt'
		rg = random_graphs.get_real_network(rand_file,None)
		if deg1_nodes == 1:
			deg1_vertices = rg.vs.select(_degree=1)
			node_indices = [v.index for v in deg1_vertices]
			rg.delete_vertices(node_indices)
		rand_PSN = build_PSN(rg)
		target = open(rand_path + network_name + '_PSN_rand_network_' + str(r) + '.txt','w')
		target.write('gene_a\tgene_b\tRAI\tJaccard_similarity\tnum_shared_partners\tnum_partners_a\tnum_partners_b\n')
		for edge in rand_PSN.es:
			gene_ids = sorted([rand_PSN.vs[edge.source]['name'],rand_PSN.vs[edge.target]['name']])
			target.write(gene_ids[0] + '\t' + gene_ids[1] + '\t' + str(edge['RAI']) + '\t' + str(edge['JS']) + '\t' + \
						 str(edge['num_shared_partners']) + '\t' + str(rg.degree(gene_ids[0])) + \
						 '\t' + str(rg.degree(gene_ids[1])) + '\n')
		target.close()
