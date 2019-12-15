# script that contains basic functions

import pandas, numpy
from custom_settings import *

# function that returns a set with all ensembl gene IDs of protein coding genes
def get_PC_genome(cursor):

	query = """select ensembl_gene_id from {}.protein_coding_genome""".format(DB_NAME)
	cursor.execute(query)
	genome = set([row[0] for row in cursor])
	return genome


# function that loads the GTEx median expression data for every PC gene and tissue restricted
# to all genes that are expressed in at least one tissue above the given expression cutoff
def get_GTEx_expr_df(connect,testis=False):

	query = """select * from {}.GTEx_expr_PC_matrix""".format(DB_NAME)
	GTEx = pandas.read_sql(query,connect)
	GTEx.set_index('ensembl_gene_id',inplace=True)
	if not testis:
		GTEx.drop('testis',axis='columns',inplace=True)
	GTEx = GTEx.loc[GTEx.max(axis=1) > GTEX_EXPR_CUTOFF,]

	return GTEx


# function that loads the TiP matrix with all values set to NaN for which the gene is not
# expressed above the given expression cutoff in the respective tissue
def get_GTEx_TiP_df(connect,testis=False):

	GTEx = get_GTEx_expr_df(connect)
	if testis:
		table_name = 'GTEx_TiP_PC_matrix'
		GTEx = get_GTEx_expr_df(connect,testis=True)
	else:
		table_name = 'GTEx_TiP_no_testis_PC_matrix'
		GTEx = get_GTEx_expr_df(connect)
	query = """select * from {}.{}""".format(DB_NAME,table_name)
	TiPmatrix = pandas.read_sql(query,connect)
	TiPmatrix.set_index('ensembl_gene_id',inplace=True)
	TiPmatrix = TiPmatrix.loc[TiPmatrix.index.isin(GTEx.index.tolist()),]
	GTEx_tissues = TiPmatrix.columns.tolist()
	for tissue in GTEx_tissues:
	    TiPmatrix.loc[TiPmatrix.index.isin(GTEx.index[GTEx[tissue]<=GTEX_EXPR_CUTOFF]),[tissue]] = numpy.NaN

	return TiPmatrix


# function that returns for every PC gene that is expressed in GTEx based on the given expression cutoff
# the max TiP value from all the tissues where the gene was observed to be expressed
def get_maxTiPvalue_series(connect,testis=False):

	TiPmatrix = get_GTEx_TiP_df(connect)
	maxTiP_series = TiPmatrix.max(axis=1,numeric_only=True)
	maxTiP_series.dropna(inplace=True)
	maxTiP_series.sort_values(inplace=True)
	return maxTiP_series
