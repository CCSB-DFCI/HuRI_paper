import pymysql
import numpy as np


#Mapping CAGE peak to entrez gene and HGNC Symbols in a dictionnary
with open ('robust_phase1_pls_2__description.txt.gz.tmp', 'r') as MappingFile:
	CAGE_peak_id_DICT = dict() #key = cagepeakID, value = [entrez-hgnc-uniprot]
	hgnc_set, uniprot_set, entrezgene_set = set(), set(), set()
	several_entrez_set = set()
    
	for i, line in enumerate(MappingFile):
		if i < 8: #HEADER of the file
			pass 

		else: #build dict: key = cage_peak_id ; value = (entrez, hgnc, uniprot)
			line = line[:-1].split('\t')

			CAGE_peak_id, entrezgene, hgnc, uniprot = line[0], line[4], line[5], line[6]
			
			if ['NA','NA','NA'] != [entrezgene, hgnc, uniprot] != ['','','']: #if a gene is related to the peak
				if 'NA' != entrezgene != '':

					if len(entrezgene.split(',')) > 1:
						#print(entrezgene)
						several_entrez_set.add(entrezgene)

					CAGE_peak_id_DICT[CAGE_peak_id] = [entrezgene, hgnc, uniprot]



#SSTAR FILE : Dictionnary to map 'ontol Id' to primary cell names
with open('Browse_samples-resource_browser.csv','r') as SstarFile:

    FFsstar_ontol_ID, Name, SampleCat, FFsstar_ontol_ID_to_PC, PC_set = set(), set(), set(), dict(), set()

    # number of line in the file: 			562 (with the header)
    # number of different primary cells: 	189

    for i, line in enumerate(SstarFile):
        line = line[1:-2].split('","')
        
        if i == 0:
            header = line
        else:
            PC_name = line[1].split(',')[0].title()
            ontol_id = line[0]
            FFsstar_ontol_ID.add(ontol_id)
            PC_set.add(PC_name)
            FFsstar_ontol_ID_to_PC[ontol_id] = PC_name



#Mapping entrez to Ensembl ids from Tables on Paros CCSB server 
#build dict (key = entrez_id, value = (enembl_id, gene type)
db = pymysql.connect('paros', user, password, 'horfeome_annotation_gencode27');
cursor = db.cursor()
cursor.execute('select distinct ensembl_gene_id_short, gene_type from gencode_transcript_annotation;')
Genetypes = cursor.fetchall()
Genetypes_dict = {line[0]:line[1] for line in Genetypes}
cursor.execute('select distinct ensembl_gene_id_short, entrez_gene_id from gencode2entrez where entrez_gene_id is not null;')
EntrezDB = cursor.fetchall()

EntrezToEnsembl_dict = {}
for i,line in enumerate(EntrezDB):
    entrez, ensembl, genetype = 'entrezgene:'+str(line[1]), line[0], Genetypes_dict[line[0]]
    if entrez not in EntrezToEnsembl_dict:
        EntrezToEnsembl_dict[entrez] = [(ensembl, genetype)]
    else:
        EntrezToEnsembl_dict[entrez].append((ensembl, genetype))



# TPM file
# Cd4s and Cd14s categories have been grouped because too similar
Cd4s = ['Cd4+ T Cells',
 'Cd4+Cd25+Cd45Ra+ Naive Regulatory T Cells',
 'Cd4+Cd25+Cd45Ra+ Naive Regulatory T Cells Expanded',
 'Cd4+Cd25+Cd45Ra- Memory Regulatory T Cells',
 'Cd4+Cd25+Cd45Ra- Memory Regulatory T Cells Expanded',
 'Cd4+Cd25-Cd45Ra+ Naive Conventional T Cells',
 'Cd4+Cd25-Cd45Ra+ Naive Conventional T Cells Expanded',
 'Cd4+Cd25-Cd45Ra- Memory Conventional T Cells',
 'Cd4+Cd25-Cd45Ra- Memory Conventional T Cells Expanded']

Cd14s = ['Cd14+ Monocyte Derived Endothelial Progenitor Cells',
 'Cd14+ Monocytes',
 'Cd14+ Monocytes - Mock Treated',
 'Cd14+ Monocytes - Treated With B-Glucan',
 'Cd14+ Monocytes - Treated With Bcg',
 'Cd14+ Monocytes - Treated With Candida',
 'Cd14+ Monocytes - Treated With Cryptococcus',
 'Cd14+ Monocytes - Treated With Group A Streptococci',
 'Cd14+ Monocytes - Treated With Ifn + N-Hexane',
 'Cd14+ Monocytes - Treated With Lipopolysaccharide',
 'Cd14+ Monocytes - Treated With Salmonella',
 'Cd14+ Monocytes - Treated With Trehalose Dimycolate (Tdm)',
 'Cd14+Cd16+ Monocytes',
 'Cd14+Cd16- Monocytes',
 'Cd14-Cd16+ Monocytes']

with open('hg19.cage_peak_phase1and2combined_tpm.osc.txt','r') as TPMFantomFile:

    URLdict = {'%20':'', '%2c':',','%3':':', '%2b':'+','%28':'(', '%29':')', '%27':'\'', '%2e':'.', '%2f':'/','%5e':''}
    mapping = dict()

    for i, line in enumerate(TPMFantomFile):
        if i > 1831:

            line = line[:-1].split('\t')
            cage_peak_id = line[0]

            if cage_peak_id in CAGE_peak_id_DICT: 
                #if cage_peak from RLE file in Sstar File : is from a sample that is a primary cell
                entrez_ids = CAGE_peak_id_DICT[cage_peak_id][0]
                for entrez_id in entrez_ids.split(','):

                    if entrez_id != '': #if there is an entrez gene id linked to the peak

                        for index in FF_ontol_indexes:
                            ontol_id = FF_ontol_ID_list[index]

                            if ontol_id in FFsstar_ontol_ID_to_PC: #is the sample linked to a primary cell?
                                pc = FFsstar_ontol_ID_to_PC[ontol_id]
                                rle = line[index]
                                

                                if entrez_id not in mapping:
                                    mapping[entrez_id] = {'Cd4+ T Cells':[], 'Cd14+ Monocytes':[]}
                                if pc not in mapping[entrez_id] and pc not in Cd4s + Cd14s:
                                    mapping[entrez_id][pc] = []
                                        
                                if pc not in Cd4s + Cd14s:
                                    mapping[entrez_id][pc].append(rle)
                                else:
                                    if pc in Cd4s:
                                        mapping[entrez_id]['Cd4+ T Cells'].append(rle)
                                    else:
                                        mapping[entrez_id]['Cd14+ Monocytes'].append(rle)


        elif i == 1831:

            for elem in URLdict:
                line = line.replace(elem, URLdict[elem])

            line = line[:-1].split('\t')
            
            FF_ontol_ID_set = {line[j][-11:] for j in range(len(line)) if j > 0} #--> size different with the list: 13 ontol id are equal
            FF_ontol_ID_list = [line[j][-11:] for j in range(len(line)) if j > 0]
            FF_ontol_set = {(line[j][-11:],line[j]) for j in range(len(line)) if j > 0} #tuple avec (ff ontol id, complete line)

            #list of all the index related to the sample (primary cell) we want to catch
            FF_ontol_indexes = []
            for j, ff_ontol_id in enumerate(FF_ontol_ID_list):
                if ff_ontol_id in FFsstar_ontol_ID: #if primary cell 
                    FF_ontol_indexes.append(j+1) #+1 because first element of line is 'gene' and not a sample
 


# write output file
def GetOrderedDictKeys(dictionary):
    '''return ordered list of dictionary keys'''
    return sorted(list(dictionary.keys()))

with open('Fantom_PC_TPM.txt', 'w') as outputfile:
    content = []
    ordered_genes_keys = GetOrderedDictKeys(mapping[list(mapping.keys())[1]])

    for i, entrez in enumerate(mapping):

        if entrez in EntrezToEnsembl_dict:

            if len(EntrezToEnsembl_dict[entrez]) == 1: 
                ensembl, genetype = EntrezToEnsembl_dict[entrez][0]  
                line = ensembl+'\t' + '\t'.join(str(np.mean(np.array(mapping[entrez][pc],dtype = float))) for pc in ordered_genes_keys)
                content.append(line)

            else: #several ensembl_ids for this entrez_id
                for (ensembl, genetype) in EntrezToEnsembl_dict[entrez]:
                    line = ensembl+'\t' + '\t'.join(str(np.mean(np.array(mapping[entrez][pc],dtype = float))) for pc in ordered_genes_keys)
                    content.append(line)


    header = 'Ensembl_IDs\t' + '\t'.join(pc for pc in ordered_genes_keys) + '\n'
    outputfile.write(header + '\n'.join(line for line in content))
