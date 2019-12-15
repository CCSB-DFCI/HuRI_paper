# script that contains code to generate functional profile matrices, some other scripts to process functional
# annotation data and some code to facilitate analysis and plotting of data by doing some pre-calculations

import numpy
import os
import sys
import pandas
import database_utils


# function that generates a matrix of all proteins in HuRI by all proteins in HuRI with their pairwise
# seq similarity
def generate_seq_similarity_matrix(cursor,path,outfile):

	print('get a mapping file between orf ids and gene ids')
	orf_gene_map_dict = {}
	query = """select distinct orf_id,ensembl_gene_id
            from horfeome_annotation_gencode27.orf_class_map_ensg
            where orf_class in ("pcORF","npcORF")"""
	cursor.execute(query)
	for row in cursor:
	    orf_id = str(row[0])
	    gene_id = row[1].split('.')[0]
	    if orf_id not in orf_gene_map_dict:
	        orf_gene_map_dict[orf_id] = []
	    orf_gene_map_dict[orf_id].append(gene_id)

	print('read in the alignment file and map to gene IDs')
	seq_align_dict = {}
	file1 = open(path + 'hi2018_paper.muscle_alignment.txt','r')
	count = 0
	for line in file1:
		tab_list = str.split(line[:-1],'\t')
		orf_id1 = tab_list[0]
		orf_id2 = tab_list[1]
		if orf_id1 != 'orf_id1':
			seq_ident = float(tab_list[2])/float(tab_list[6])
			if orf_id1 in orf_gene_map_dict and orf_id2 in orf_gene_map_dict:
				for gene_a in orf_gene_map_dict[orf_id1]:
					for gene_b in orf_gene_map_dict[orf_id2]:
						pair = tuple(sorted([gene_a,gene_b]))
						if pair not in seq_align_dict:
							seq_align_dict[pair] = 0
						seq_align_dict[pair] = max(seq_align_dict[pair],seq_ident)
		count += 1
		if count % 100000 == 0:
			print(count)
	file1.close()

	print('get all gene ids')
	gene_ids = sorted(list(set([tup[0] for tup in seq_align_dict.keys()] + [tup[1] for tup in seq_align_dict.keys()])))

	print('create a matrix of seq similarities')
	matrix = numpy.ones((len(gene_ids),len(gene_ids)),dtype=float)
	for pair,seq_ident in seq_align_dict.items():
		index_a = gene_ids.index(pair[0])
		index_b = gene_ids.index(pair[1])
		matrix[index_a,index_b] = seq_ident
		matrix[index_b,index_a] = seq_ident

	print('write out the matrix')
	target = open(outfile,'w')
	target.write('\t' + '\t'.join(gene_ids) + '\n')
	for i in range(len(gene_ids)):
		target.write(gene_ids[i] + '\t' + '\t'.join([str(v) for v in matrix[i,]]) + '\n')
	target.close()


# function that processes the co-fitness profile matrix by mapping it to Ensembl gene IDs and dealing with
# changes to the dataset because of the mapping
def map_co_fitness_matrix_to_gene_ids(cursor,co_fitness_matrix_file,outfile):

	# get a mapping dict between gene symbols and gene IDs
	query = """select distinct ensembl_gene_id_short,symbol from horfeome_annotation_gencode27.gencode2entrez"""
	cursor.execute(query)
	symbol_id_map_dict = {}
	for row in cursor:
		gene_id = row[0]
		symbol = row[1]
		symbol_id_map_dict[symbol] = gene_id

	file1 = open(co_fitness_matrix_file,'r')
	for line in file1:
		gene_names = line[:-1].split('\t')[1:]
		break
	matrix = numpy.loadtxt(co_fitness_matrix_file,dtype='float',delimiter='\t',skiprows=1,usecols=range(1,len(gene_names)+1))

	gene_ids = []
	indices = []
	for i,gene_name in enumerate(gene_names):
		if gene_name in symbol_id_map_dict:
			gene_ids.append(symbol_id_map_dict[gene_name])
			indices.append(i)
		else:
			gene_ids.append(None)

	target = open(outfile,'w')
	target.write('\t' + '\t'.join(filter(lambda v: v is not None,gene_ids)) + '\n')
	for i in indices:
		target.write(gene_ids[i] + '\t' + '\t'.join([str(matrix[i,v]) for v in indices]) + '\n')
	target.close()


# function that reads in the files from SEEK DB to generate one PCC matrix restricted to proteins in HI-union
def generate_SEEK_matrix(SEEK_path,cursor,outfile_prefix):

	# get the set of genes to be included in the expression matrix
	query = """select distinct a.ensembl_gene_id from
	horfeome_annotation_gencode27.orf_class_map_ensg a,
	horfeome_annotation_gencode27.gencode_transcript_annotation b,
	(select distinct ad_orf_id orf_id from hi_ref_retest.retest where standard_batch in ("Hs15","Hs14") and final_score="1"
	union
	(select distinct db_orf_id from hi_ref_retest.retest where standard_batch in ("Hs15","Hs14") and final_score="1")) as c
	where a.in_new_space_3=1 and a.orf_class='pcORF' and a.ensembl_gene_id=b.ensembl_gene_id and b.gene_type="protein_coding" and a.orf_id=c.orf_id
	union
	(select distinct a.ensembl_gene_id from
	(select interactor_a ensembl_gene_id from
	hi_ref.huri_exp_info
	union
	(select interactor_b from
	hi_ref.huri_exp_info
	)) as a,
	horfeome_annotation_gencode27.gencode_transcript_annotation b,
	horfeome_annotation_gencode27.orf_class_map_ensg c
	where a.ensembl_gene_id=b.ensembl_gene_id and b.gene_type='protein_coding' and a.ensembl_gene_id=c.ensembl_gene_id and c.orf_class='pcORF')"""
	cursor.execute(query)
	huri_genes = set([row[0].split('.')[0] for row in cursor])
	print('Number of HI-union genes:', len(huri_genes))

	# get a mapping between gene IDs and Ensembl gene ids
	query = """select distinct entrez_gene_id,ensembl_gene_id_short from horfeome_annotation_gencode27.gencode2entrez"""
	cursor.execute(query)
	map_dict = {}
	for row in cursor:
		map_dict[str(row[0])] = row[1]
	print('Number of geneIDs mapped to Ensg IDs:', len(map_dict))

	# get set of gene IDs used in SEEK as intersection between file names and row names in one file and map those to ensemble gene ids
	files = os.listdir(SEEK_path)
	gene_ids = set([s.split('.')[0] for s in files])
	file1 = open(SEEK_path + files[0],'r')
	entries = file1.readlines()
	gene_ids = list(gene_ids.intersection(set([line.split('\t')[0] for line in entries])))
	matched_gene_ids = []
	matched_ensg_ids = []
	for gene_id in gene_ids:
		if gene_id in map_dict and map_dict[gene_id] in huri_genes and map_dict[gene_id] not in matched_ensg_ids:
			matched_gene_ids.append(gene_id)
			matched_ensg_ids.append(map_dict[gene_id])
	print('Number of retained ensg IDs:', len(matched_ensg_ids))

	# get a dict that maps from gene ID to index in matrix
	id_index_dict = dict([(matched_gene_ids[i],i) for i in range(len(matched_gene_ids))])

	# fill a matrix with the PCCs for the selected list of gene IDs
	print('Start filling the PCC matrix')
	matrix = numpy.zeros((len(matched_ensg_ids),len(matched_ensg_ids)),dtype=float)
	for f,file_name in enumerate(files):
		if f % 1000 == 0:
			print(f)
		gene_a = file_name.split('.')[0]
		if gene_a in id_index_dict:
			file1 = open(SEEK_path + file_name,'r')
			entries = file1.readlines()
			file1.close()
			for line in entries:
				tab_list = str.split(line[:-1],'\t')
				gene_b = tab_list[0]
				pcc = float(tab_list[1])
				if gene_b in id_index_dict:
					matrix[id_index_dict[gene_a],id_index_dict[gene_b]] = pcc
					matrix[id_index_dict[gene_b],id_index_dict[gene_a]] = pcc

	# write out
	print('Write text')
	outfile_text = outfile_prefix + '.txt'
	target = open(outfile_text,'w')
	target.write('\t' + '\t'.join(matched_ensg_ids) + '\n')
	for i in range(len(matched_ensg_ids)):
		target.write(matched_ensg_ids[i] + '\t' + '\t'.join([str(v) for v in matrix[i,:]]) + '\n')
	target.close()
	print('Write binary')
	outfile_binary = outfile_prefix + '.npy'
	numpy.save(outfile_binary,matrix)


def count_edges_per_cutoff(nw_name,num_rand,inpath,outpath):

	cutoffs = [0.001,0.01] + [i/10.0 for i in range(1,11)]

	count_matrix = numpy.zeros((num_rand,len(cutoffs)),dtype='int')
	for r in range(num_rand):
	    if r % 10 == 0:
	        print(r)
	    network_file = inpath + nw_name + '/' + nw_name + '_rand_network_' + str(r) + '.txt'
	    PSN = pandas.read_table(network_file,header=0)
	    for c,cutoff in enumerate(cutoffs):
	        count_matrix[r][c] = PSN.loc[PSN['Jaccard_similarity'] >= cutoff,].shape[0]

	outfile = outpath + nw_name + '_rand_edge_count_per_JScutoff.npy'
	numpy.save(outfile,count_matrix)


if __name__ == '__main__':

	connect = database_utils.get_connection()
	cursor = connect.cursor()
	path = '../data/katjas_data/PSN_analysis/'

	mode = sys.argv[1]

	if mode == '1':
		outfile = path + 'seq_similarity_matrix.txt'
		generate_seq_similarity_matrix(cursor,path,outfile)
	elif mode == '2':
		co_fitness_matrix_file = path + 'avana_2017_dep_corr.tsv'
		outfile = path + 'avana_2017_dep_corr_ENSG_ID.tsv'
		map_co_fitness_matrix_to_gene_ids(cursor,co_fitness_matrix_file,outfile)
	elif mode == '4':
		SEEK_path = path + 'results.human.single.gene.cor/'
		outfile_prefix = path + 'SEEK_matrix'
		generate_SEEK_matrix(SEEK_path,cursor,outfile_prefix)
	elif mode == '6':
		nw_name = sys.argv[2]
		num_rand = int(sys.argv[3])
		inpath = path
		outpath = path
		count_edges_per_cutoff(nw_name,num_rand,inpath,outpath)
