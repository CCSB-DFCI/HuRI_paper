# script that contains objects to load and query disease and edgotyping data

import sys
import database_utils
import edgotype_classes


class Gene(edgotype_classes.Gene):

	def __init__(self,gene_id):

		super(Gene, self).__init__(gene_id)
		self.orf_ids = []
		self.PPIs = {}

class Allele(edgotype_classes.Allele):

	def __init__(self,allele_id):

		super(Allele, self).__init__(allele_id)

class HI_Interaction():

	def __init__(self,gene_a,gene_b):

		self.gene_a = gene_a
		self.gene_b = gene_b
		self.orf_id_a = ''
		self.orf_id_b = ''
		self.version1 = 0
		self.version2 = 0
		self.version3 = 0
		self.v1_ori = ''
		self.v2_ori = ''
		self.v3_ori = ''

class Allele_Interaction(edgotype_classes.Interaction):

	def __init__(self,allele_id,partner_id,assay_id):

		super(Allele_Interaction, self).__init__(allele_id,partner_id,assay_id)

class Disease(edgotype_classes.Disease):

	def __init__(self,omim_id):

		super(Disease, self).__init__(omim_id)
		self.omim_id = omim_id
		self.name = ''
		self.ps_info = set()
		self.tissues = set()
		self.causal_genes = set()
		self.CG_tissue_pairs = set()

class Analysis(edgotype_classes.Analysis):

	def __init__(self,name):

		super(Analysis, self).__init__(name)
		self.name = name
		self.genes = {}
		self.alleles = {}
		self.hi_interactions = {}
		self.allele_interactions = {}
		self.diseases = {}
		self.tissues = []

	# function that loads the gene information into the analysis object
	def load_genes(self,cursor,tissues):

		# get all genes that are in HI-III
		query = """select distinct ensembl_gene_id_a from hi2018_paper.hi_iii union (select distinct ensembl_gene_id_b from hi2018_paper.hi_iii)"""
		cursor.execute(query)
		for row in cursor:
			gene_id = row[0]
			gene_obj = Gene(gene_id)
			self.genes[gene_id] = gene_obj

		# add their ORF IDs
		query = """select distinct orf_id,ensembl_gene_id
				from hi2018_paper.orf_class_map_ensg
				where in_new_space_3=1 and orf_class='pcORF'"""
		cursor.execute(query)
		for row in cursor:
			gene_id = row[1].split('.')[0]
			orf_id = int(row[0])
			if gene_id in self.genes:
				self.genes[gene_id].orf_ids.append(orf_id)

		# get the symbol and description for these genes
		query = """select distinct a.ensembl_gene_id,a.symbol,b.DESCRIPTION
				from hi2018_paper.orf_class_map_ensg a,
				entrez_gene.gene_info_human_20170914 b
				where a.entrez_gene_id=b.GENE_ID"""
		cursor.execute(query)
		for row in cursor:
			gene_id = row[0].split('.')[0]
			if gene_id in self.genes:
				self.genes[gene_id].symbol = row[1]
				self.genes[gene_id].descr = row[2]

		# get the expr and TiP values for these genes
		query_tissues = ['`' + t + '`' for t in tissues]
		cols = ','.join(['ensembl_gene_id'] + query_tissues)
		print(cols)
		query = """select {} from hi2018_paper.GTEx_expr_PC_matrix""".format(cols)
		cursor.execute(query)
		for row in cursor:
			gene_id = row[0]
			if gene_id in self.genes:
				self.genes[gene_id].tissue_expr_values = [float(v) for v in row[1:]]
		query = """select {} from hi2018_paper.GTEx_TiP_no_testis_PC_matrix""".format(cols)
		cursor.execute(query)
		for row in cursor:
			gene_id = row[0]
			if gene_id in self.genes:
				self.genes[gene_id].tissue_TiP_values = [float(v) for v in row[1:]]

		# get the diseases
		if len(self.diseases) == 0:
			print('Disease data not loaded.')
		else:
			for omim_id,disease_obj in self.diseases.items():
				for gene_id in disease_obj.causal_genes:
					if gene_id in self.genes:
						self.genes[gene_id].diseases[omim_id] = disease_obj


	# get all disease alleles loaded from the experiment and annotate them with more information
	def load_disease_alleles(self,cursor):

		std_batches_str = ','.join(["'" + v + "'" for v in self.std_batch_dict.keys()])

		# add all disease alleles
		query = """select distinct b.db_mut_id,min(c.ensembl_gene_id)
				from cegs_retest.retest b,
				(select distinct orf_id,ensembl_gene_id
					from hi2018_paper.orf_class_map_ensg
					where orf_class='pcORF' and in_new_space_3=1) as c
				where b.standard_batch in ({}) and b.db_mut_id!=0
				and b.db_orf_id=c.orf_id group by b.db_mut_id""".format(std_batches_str)
		cursor.execute(query)
		for row in cursor:
			allele_id = int(row[0])
			gene_id = row[1].split('.')[0]
			if gene_id in self.genes:
				allele_obj = Allele(allele_id)
				allele_obj.gene_id = gene_id
				self.alleles[allele_id] = allele_obj

		# add to alleles information on amino acid change and associated tissue-specific diseases
		query = """select distinct mut_id,hgvs,PhenotypeIDS,ClinicalSignificance
					from cegs.CEGS_ClinVar_2018"""
		cursor.execute(query)
		for row in cursor:
			allele_id = int(row[0])
			if allele_id in self.alleles:
				self.alleles[allele_id].mutation = row[1]
				annot = row[2]
				annot_words = annot.split(';')
				omim_ids = set()
				for annot_word in annot_words:
					id_words = annot_word.split(',')
					for id_word in id_words:
						if id_word.find('OMIM:') > -1:
							omim_ids.add(int(id_word.split(':')[1]))
				for omim_id in omim_ids:
					if omim_id in self.diseases:
						self.alleles[allele_id].diseases[omim_id] = self.diseases[omim_id]
				mut_class = row[3]
				self.alleles[allele_id].allele_type = mut_class

		# link disease alleles to gene objects
		for allele_id,allele_obj in self.alleles.items():
			self.genes[allele_obj.gene_id].disease_alleles[allele_id] = allele_obj


	# function that loads the HI-III data
	def load_hi_interactions(self,cursor):

		query = """select ensembl_gene_id_a,ensembl_gene_id_b,orientation,assay_id,orf_id_a,orf_id_b
				from hi2018_paper.huri_exp_info where dataset='validated'
				"""
		cursor.execute(query)
		for row in cursor:
			gene_a = row[0]
			gene_b = row[1]
			ori = row[2]
			assay_id = int(row[3])
			orf_id_a = int(row[4])
			orf_id_b = int(row[5])
			if (gene_a,gene_b) not in self.hi_interactions:
				int_obj = HI_Interaction(gene_a,gene_b)
				int_obj.orf_id_a = orf_id_a
				int_obj.orf_id_b = orf_id_b
				self.hi_interactions[(gene_a,gene_b)] = int_obj
				self.genes[gene_a].PPIs[gene_b] = int_obj
				self.genes[gene_b].PPIs[gene_a] = int_obj
			int_obj = self.hi_interactions[(gene_a,gene_b)]
			if assay_id == 1:
				int_obj.version1 = 1
				int_obj.v1_ori = ori
			if assay_id == 2:
				int_obj.version2 = 1
				int_obj.v2_ori = ori
			if assay_id == 6:
				int_obj.version3 = 1
				int_obj.v3_ori = ori


	# function that loads and processes the edgotyping data
	def load_allele_interactions(self,cursor):

		std_batches_str = ','.join(["'" + v + "'" for v in self.std_batch_dict.keys()])

		# for every causal gene get a list of all valid interaction partners that worked in this assay with
		# the corresponding growth score
		# get a temporary dictionary that maps between orf_ids and gene_ids
		temp_dict = {}
		for gene_id,gene_obj in self.genes.items():
			for orf_id in gene_obj.orf_ids:
				temp_dict[orf_id] = gene_id
		query = """select distinct a.ad_orf_id,a.db_orf_id,c.growth_test,a.standard_batch
				from cegs_retest.retest a,
				cegs_retest.retest_scoring c,
				cegs_retest.retest_source d
				where a.standard_batch in ({}) and a.standard_batch=d.standard_batch and a.ad_orf_id=d.ad_orf_id and
				a.db_orf_id=d.db_orf_id and a.db_mut_id=d.db_mut_id and a.ad_mut_id=d.ad_mut_id and d.source!='RRS'
				and d.source!='rrs' and a.db_mut_id=0 and a.final_score='1' and a.final_score_id=c.score_id and
				c.growth_test>1 and a.seq_confirmation_db='y' and a.seq_confirmation_lw='y'""".format(std_batches_str)
		cursor.execute(query)
		rows = cursor.fetchall()
		for row in rows:
			ad_orf_id = int(row[0])
			db_orf_id = int(row[1])
			gs = int(row[2])
			assay_id = self.std_batch_dict[row[3]]
			if ad_orf_id in temp_dict and db_orf_id in temp_dict:
				ad_gene_id = temp_dict[ad_orf_id]
				db_gene_id = temp_dict[db_orf_id]
				self.genes[db_gene_id].wt_PPIs[(ad_gene_id,assay_id)] = gs

		# get the interaction data for the mutants
		query = """select distinct a.ad_orf_id,a.db_orf_id,a.db_mut_id,c.growth_test,a.final_score,a.standard_batch
				from cegs_retest.retest a,
				cegs_retest.retest_scoring c,
				cegs_retest.retest_source d
				where a.standard_batch in ({}) and a.standard_batch=d.standard_batch and a.ad_orf_id=d.ad_orf_id and
				a.db_orf_id=d.db_orf_id and a.db_mut_id=d.db_mut_id and a.ad_mut_id=d.ad_mut_id and d.source!='RRS'
				and d.source!='rrs' and a.db_mut_id!=0 and a.final_score_id=c.score_id and a.seq_confirmation_db='y' and
				a.seq_confirmation_lw='y'""".format(std_batches_str)
		cursor.execute(query)
		rows = cursor.fetchall()
		# for every line in the retest table
		for row in rows:
			# if this corresponds to a causal gene in my list and the partner is valid
			ad_orf_id = int(row[0])
			db_orf_id = int(row[1])
			allele_id = int(row[2])
			score = row[4]
			assay_id = self.std_batch_dict[row[5]]
			if ad_orf_id in temp_dict and db_orf_id in temp_dict:
				ad_gene_id = temp_dict[ad_orf_id]
				db_gene_id = temp_dict[db_orf_id]
				if (ad_gene_id,assay_id) in self.genes[db_gene_id].wt_PPIs:
					# if the interaction worked and there is a valid growth score
					if (score == '1' or score == '0') and row[3] is not None and row[3] is not 'NA':
						# get the growth score from wt and mut allele, subtract, and determine
						# whether the edge is perturbed
						gs_wt = self.genes[db_gene_id].wt_PPIs[(ad_gene_id,assay_id)]
						gs_mut = int(row[3])
						gs_diff = gs_wt - gs_mut
						if gs_diff >= 2:
							perturbed = 1
						else:
							perturbed = 0

						# make an interaction object, fill it and link it with the allele object
						int_obj = Allele_Interaction(allele_id,ad_gene_id,assay_id)
						int_obj.score = score
						int_obj.GS_decrease = gs_diff
						int_obj.GS_wt = gs_wt
						int_obj.GS_mut = gs_mut
						int_obj.perturbed = perturbed
						self.alleles[allele_id].interactions[(allele_id,ad_gene_id,assay_id)] = int_obj
						self.allele_interactions[(allele_id,ad_gene_id,assay_id)] = int_obj



def get_hi_iii_CG_TiP_data(std_batches):

	connect = database_utils.get_connection()
	cursor = connect.cursor()

	disease_table = 'hi2018_paper.OMIM_disease_tissue'

	file1 = open('../data/katjas_data/GTEx/expression_datasets/GTExv6_JPnorm_subtypes_no_cell_no_testis.txt','r')
	tissues = [line[:-1] for line in file1.readlines()]
	file1.close()

	analysis1 = Analysis('analysis1')
	analysis1.tissues = tissues
	analysis1.get_std_batch_assay_id_dict(cursor,std_batches)
	analysis1.load_diseases(cursor,disease_table)
	print('Number of diseases loaded:',len(analysis1.diseases))

	analysis1.load_genes(cursor,tissues)
	print('Number of genes loaded:',len(analysis1.genes))

	analysis1.load_disease_alleles(cursor)
	print('Number of alleles loaded:',len(analysis1.alleles))

	analysis1.load_hi_interactions(cursor)
	print('Number of HI-III interactions loaded:',len(analysis1.hi_interactions))

	analysis1.load_allele_interactions(cursor)
	print('Number of allele interactions loaded:',len(analysis1.allele_interactions))

	return analysis1
