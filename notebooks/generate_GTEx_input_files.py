# script to generate the GTEx expression and TiP files that are used as input for other scripts

import os
import pandas, numpy
import database_utils, basic_functions
from custom_settings import *

def generate_tissue_name_table(cursor,testis=False):

	files = os.listdir('../external_data/GTEx_sample_data/')
	tissues = []
	for file in files:
		tissue = file.split('.')[0]
		if tissue[:5] != 'cells' and (tissue != 'testis' or (testis and tissue == 'testis')):
			tissues.append(tissue)
	tissues.sort()

	query = """delete from {}.GTEx_tissue_names""".format(DB_NAME)
	cursor.execute(query)

	for tissue in tissues:
		query = """insert into {}.GTEx_tissue_names (tissue_names) values (%s)""".format(DB_NAME)
		cursor.execute(query,(tissue,))


def generate_GTEx_expression_table(PC_genome,cursor):

	files = os.listdir('../external_data/GTEx_sample_data/')
	files.sort()

	# initialize GTEx DataFrame with median expression of all genes in first tissue
	tissue = files[0].split('.')[0]
	df = pandas.read_table('../external_data/GTEx_sample_data/' + files[0],index_col=0)
	gtex = pandas.DataFrame({tissue:df.median(axis=1)},index=df.index)
	gtex.index.rename('ensembl_gene_id',inplace=True)

	# add to GTEx DataFrame median expression of all genes in all remaining tissues
	for file in files[1:]:
		tissue = file.split('.')[0]
		if tissue[:5] != 'cells':
			df = pandas.read_table('../external_data/GTEx_sample_data/' + file,index_col=0)
			gtex[tissue] = df.median(axis=1)

	# remove ensembl gene ID version information and any non-coding genes from GTEx DataFrame
	gtex.rename(mapper=(lambda x: x.split('.')[0]),axis='index',inplace=True)
	gtex = gtex.loc[gtex.index.isin(PC_genome),]

	# upload GTEx DataFrame to DB table
	query = """delete from {}.GTEx_expr_PC_matrix""".format(DB_NAME)
	cursor.execute(query)
	for i in range(gtex.shape[0]):
	    query = 'insert into {}.GTEx_expr_PC_matrix values ('.format(DB_NAME) + ','.join(['%s' for c in range(gtex.shape[1]+1)]) + ')'
	    cursor.execute(query,tuple([gtex.index[i]] + gtex.iloc[i,].tolist()))


def generate_TiP_value_table(PC_genome,connect,testis=False):

	# read in all the sample expression data into one DataFrame
	files = os.listdir('../external_data/GTEx_sample_data/')
	files.sort()

	# initialize sample DataFrame with samples from first tissue
	tissue = files[0].split('.')[0]
	gtex = pandas.read_table('../external_data/GTEx_sample_data/' + files[0],index_col=0)
	gtex.index.rename('ensembl_gene_id',inplace=True)

	# add to GTEx DataFrame samples of all remaining tissues
	for file in files[1:]:
		tissue = file.split('.')[0]
		if tissue[:5] != 'cells' and (tissue != 'testis' or (testis and tissue == 'testis')):
			df = pandas.read_table('../external_data/GTEx_sample_data/' + file,index_col=0)
			df.index.rename('ensembl_gene_id',inplace=True)
			gtex = gtex.merge(df,left_index=True,right_index=True)

	# remove ensembl gene ID version information and any non-coding genes from samples DataFrame
	gtex.rename(mapper=(lambda x: x.split('.')[0]),axis='index',inplace=True)
	gtex = gtex.loc[gtex.index.isin(PC_genome),]

	stats_df = pandas.DataFrame({'all_median':gtex.median(axis=1),
								 'all_IQR':gtex.quantile(q=0.75,axis=1) - gtex.quantile(q=0.25,axis=1)},
								 index=gtex.index)

	query = """select * from {}.GTEx_expr_PC_matrix""".format(DB_NAME)
	expr_df = pandas.read_sql(query,connect)
	expr_df.set_index('ensembl_gene_id',inplace=True)
	if not testis:
		expr_df.drop('testis',axis='columns',inplace=True)

	TiP_df = pandas.DataFrame({},index=expr_df.index)
	col_names = sorted(expr_df.columns.tolist())
	for tissue in col_names:
		TiP_df[tissue] = (expr_df[tissue] - stats_df['all_median'])/stats_df['all_IQR']

	# upload TiP DataFrame to DB table
	cursor = connect.cursor()
	if not testis:
		table_name = 'GTEx_TiP_no_testis_PC_matrix'
	else:
		table_name = 'GTEx_TiP_PC_matrix'

	query = """delete from {}.{}""".format(DB_NAME,table_name)
	cursor.execute(query)

	for i in range(TiP_df.shape[0]):
	    query = 'insert into {}.{} values ('.format(DB_NAME,table_name) + ','.join(['%s' for c in range(TiP_df.shape[1]+1)]) + ')'
	    cursor.execute(query,tuple([TiP_df.index[i]] + TiP_df.iloc[i,].tolist()))


def get_maxTiPvalue_table(connect):

	TiPmatrix = basic_functions.get_GTEx_TiP_df(connect)
	print(TiPmatrix.shape)
	cursor = connect.cursor()
	query = """delete from {}.GTEx_maxTiP_info""".format(DB_NAME)
	cursor.execute(query)

	df = pandas.DataFrame({'maxTiP':TiPmatrix.max(axis=1),'maxTissue':TiPmatrix.idxmax(axis=1)},index=TiPmatrix.index)
	print(df.shape)

	for i in range(df.shape[0]):
		query = """insert into {}.GTEx_maxTiP_info values(%s,%s,%s)""".format(DB_NAME)
		cursor.execute(query,(df.index[i],df.ix[i,'maxTiP'],df.ix[i,'maxTissue']))


if __name__ == '__main__':

	connect = database_utils.get_connection()
	cursor = connect.cursor()

	# generate DB table with tissue names
	print('Generate tissue name table')
	generate_tissue_name_table(cursor)

	# generate GTEx expression DB table
	PC_genome = basic_functions.get_PC_genome(cursor)
	print('Generate GTEx expression table')
	generate_GTEx_expression_table(PC_genome,cursor)

	# generate TiP value DB table
	print('Generate GTEx TiP table')
	generate_TiP_value_table(PC_genome,connect)
	generate_TiP_value_table(PC_genome,connect,testis=True)

	# generate max TiP value DB table
	get_maxTiPvalue_table(connect)
