# function that generates HI-III used for the PSN analyses

import huri_db_utils

if __name__ == '__main__':

	path = '../data/katjas_data/PSN_analysis/'

	# generate a file for HI-III, read in with homodimers but write out without homodimer ppis
	hi_iii = huri_db_utils.load_nw_hi_iii(id_type='ensembl_gene_id',fmt='igraph',cache='off')
	target = open(path + 'HI-III.txt','w')
	target.write('gene_a\tgene_b\n')
	for edge in hi_iii.es:
		gene_a = hi_iii.vs[edge.source]['name']
		gene_b = hi_iii.vs[edge.target]['name']
		if gene_a != gene_b:
			if gene_a < gene_b:
				target.write(gene_a + '\t' + gene_b + '\n')
			else:
				target.write(gene_b + '\t' + gene_a + '\n')
	target.close()
