# script that contains objects to load and query disease and edgotyping data

import sys
import database_utils

class Gene(object):

	def __init__(self,gene_id):

		self.gene_id = gene_id
		self.symbol = ''
		self.descr = ''
		self.orf_id = ''
		self.disease_alleles = {}
		self.tissue_expr_values = []
		self.tissue_TiP_values = []
		self.diseases = {}
		self.wt_PPIs = {}

class Allele(object):

	def __init__(self,allele_id):

		self.allele_id = allele_id
		self.allele_type = None
		self.gene_id = ''
		self.mutation = ''
		self.diseases = {}
		self.interactions = {}

class Interaction(object):

	def __init__(self,allele_id,partner_id,assay_id):

		self.allele_id = allele_id
		self.partner_id = partner_id
		self.assay_id = assay_id
		self.score = None
		self.perturbed = None
		self.GS_decrease = None
		self.GS_wt = None
		self.GS_mut = None

class Disease(object):

	def __init__(self,omim_id):

		self.omim_id = omim_id
		self.name = ''
		self.ps_info = set()
		self.tissues = set()
		self.causal_genes = set()
		self.CG_tissue_pairs = set()


class Analysis(object):

	def __init__(self,name):

		self.name = name
		self.genes = {}
		self.alleles = {}
		self.interactions = {}
		self.diseases = {}
		self.tissues = []
		self.std_batch_dict = {}

	# function that loads disease information into the analysis object
	def load_diseases(self,cursor,db_table):

		query = """select distinct ensembl_gene_id,disease_tissue,disease,omim_id,phenotypic_series,ps_id
				from {}
				where number_tissues<=3 and gene_type='protein_coding' and
					disease_tissue!='testis' and tissue_expr_value>5""".format(db_table)
		cursor.execute(query)

		for row in cursor:
			omim_id = int(row[3])
			if omim_id not in self.diseases:
				disease_obj = Disease(omim_id)
				disease_obj.name = row[2]
				self.diseases[omim_id] = disease_obj
			disease_obj = self.diseases[omim_id]
			disease_obj.causal_genes.add(row[0])
			disease_obj.tissues.add(row[1])
			disease_obj.ps_info.add((row[5],row[4]))
			disease_obj.CG_tissue_pairs.add((row[0],row[1]))

	# function that generates and returns a dict: std_batch -> assay_id
	def get_std_batch_assay_id_dict(self,cursor,std_batches):

		std_batches_str = ','.join(["'" + v + "'" for v in std_batches])
		query = """select distinct standard_batch,assay_id
					from cegs_retest.retest_batch
					where standard_batch in ({})""".format(std_batches_str)
		cursor.execute(query)
		self.std_batch_dict = dict([(row[0],int(row[1])) for row in cursor.fetchall()])
