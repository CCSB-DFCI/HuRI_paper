# script to prepare the functional co-annotation files

import pandas
from huri_external_data import *
from huri_db_utils import *
import database_utils

def format_reactome(outfile):

	df_pw = load_reactome_gene_to_pathway()
	pw_genes = df_pw['ensembl_gene_id'].unique()
	annot_dict = {}
	for gene_id in pw_genes:
		annot_dict[gene_id] = set(df_pw.loc[df_pw['ensembl_gene_id'] == gene_id,'reactome_pathway_id'].values)

	target = open(outfile,'w')
	target.write('ensembl_gene_id\tpathway_ids\n')
	for gene_id,pathway_ids in annot_dict.items():
		target.write(gene_id + '\t' + '|'.join(sorted(list(pathway_ids))) + '\n')
	target.close()


def format_cellular_localization_data(outfile):

	df_cl = subcellular_location()
	gene_ids = df_cl.index.values
	cl_names = df_cl.columns.values
	annot_dict = dict([(gene_id,set()) for gene_id in gene_ids])
	for cl_name in cl_names:
		genes = df_cl.loc[df_cl[cl_name],].index.values
		for gene in genes:
			annot_dict[gene].add(cl_name[3:].replace(' ','_'))

	target = open(outfile,'w')
	target.write('ensembl_gene_id\tcompartments\n')
	for gene_id,cl_names in annot_dict.items():
		target.write(gene_id + '\t' + '|'.join(sorted(list(cl_names))) + '\n')
	target.close()


def format_bioplex_complexes(infile,outfile):

	connect = database_utils.get_connection()
	df_bp = pandas.read_table(infile)
	query = """select distinct ensembl_gene_id_short,entrez_gene_id
	           from horfeome_annotation_gencode27.gencode2entrez where entrez_gene_id is not NULL"""
	geneid_ensg_map = pandas.read_sql(query,connect)
	df_bp = df_bp.merge(geneid_ensg_map,left_on='GeneID',right_on='entrez_gene_id',how='inner').drop(['entrez_gene_id'],axis=1)
	gene_ids = df_bp['ensembl_gene_id_short'].unique()
	annot_dict = {}
	for gene_id in gene_ids:
		annot_dict[gene_id] = set(df_bp.loc[df_bp['ensembl_gene_id_short'] == gene_id,'Cluster Number'].values)
		if len(annot_dict[gene_id]) > 1:
			print gene_id

	target = open(outfile,'w')
	target.write('ensembl_gene_id\tBioplex_complex_ids\n')
	for gene_id,cluster_ids in annot_dict.items():
		target.write(gene_id + '\t' + '|'.join([str(i) for i in sorted(list(cluster_ids))]) + '\n')
	target.close()


if __name__ == '__main__':

	path = '../data/katjas_data/PSN_analysis/'

	outfile = path + 'fct_annotations/reactome_pathway_info.txt'
	format_reactome(outfile)

	outfile = path + 'fct_annotations/cell_atlas_localization_info.txt'
	format_cellular_localization_data(outfile)

	infile = path + 'fct_annotations/Bioplex2.0_cluster_suppl_table7.txt'
	outfile = path + 'fct_annotations/Bioplex2.0_complex_info.txt'
	format_bioplex_complexes(infile,outfile)
