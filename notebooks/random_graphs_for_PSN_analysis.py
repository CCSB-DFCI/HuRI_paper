
import sys, os
import random_graphs

if __name__ == '__main__':

	path = sys.argv[1]
	network_name = sys.argv[2]
	num_graphs = int(sys.argv[3])

	outpath = path + network_name + '/'
	if not os.path.exists(outpath):
		os.mkdir(outpath)

	network_file = path + network_name + '.txt'
	random_graphs.get_random_graphs(network_file,outpath,network_name,num_graphs,None)
