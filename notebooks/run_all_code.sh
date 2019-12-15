# executable shell script to document and run the various scripts

"""
code to get the GTEx data
"""

# download and save R object with GTEx v6 data, processed and normalized as described in Paulson et al. BMC Bioinformatics 2017
Rscript download_GTEx_data.R
# generate a file for every tissue containing the expression data for every sample that was sequenced
Rscript GTEx_generate_file_per_subtype.R
# rename some tissues and corresponding files
mv ../external_data/adipose_visceral_\(omentum\).txt ../external_data/adipose_visceral_omentum.txt
mv ../external_data/brain-0.txt ../external_data/brain_basal_ganglia.txt
mv ../external_data/brain-1.txt ../external_data/brain_cerebellum.txt
mv ../external_data/brain-2.txt ../external_data/brain_other.txt

# generate DB tables for tissues to be used, median expression of all PC genes in all tissues, and TiP values
python3 generate_GTEx_input_files.py

# generate plots that characterize TiP genes and impact on set of TiP genes when removing testis samples prior to
# TiP gene calculation (Ext. Data Fig 10)
jupyter nbconvert --to notebook --execute GTEx_TiP_characterization.ipynb
rm GTEx_TiP_characterization.nbconvert.ipynb


"""
code to process FANTOM and HPA data and make space plots in Fig 1
"""
# generate plots that characterize TiP genes and impact on set of TiP genes when removing testis samples prior to
# TiP gene calculation (Ext. Data Fig 10)
jupyter nbconvert --to notebook --execute space_descr_plots.ipynb
rm space_descr_plots.nbconvert.ipynb


"""
code to do all analyses with CORUM (Ext Data Fig. 2i and 3)
"""
jupyter nbconvert --to notebook --execute corum_huri.ipynb
rm corum_huri.nbconvert.ipynb


"""
code to generate PSN and random PSNs of HuRI and to perform analyses on PSN
"""
# get HuRI without homodimeric PPIs
python3 ./network_generation.py
# randomize that network
python3 ./random_graphs_for_PSN_analysis.py ../data/katjas_data/PSN_analysis/ HI-III 1000
# get PSN and random PSNs from HuRI
python3 ./generate_PSNs.py HI-III True 0 99 ../data/katjas_data/PSN_analysis/ 1
# compute number of edges in PSN for increasing Jaccard similarity cutoffs
python3 ./fct_network_functions.py 6 HI-III_PSN 100
# process gene function annotation files for pathway, protein complex, and subcellular compartment data
python3 ./fct_annotation_files.py
# compute functional enrichments for pathway, protein complex, and subcellular compartment data
python3 ./PPI_coannot_overlap_calculations.py HI-III 100 PPI
python3 ./PPI_coannot_overlap_calculations.py HI-III_PSN 100 Jaccard_similarity
# compute functional PSNs
python3 ./fct_network_functions.py 1
python3 ./fct_network_functions.py 2
python3 ./fct_network_functions.py 4
# compute functional enrichments for seq identity, co-fitness and co-expression data
python3 ./PPI_fctPSN_overlap_calculations.py HI-III avana_2017_dep_corr_ENSG_ID coFitnessPSN 100 PPI 0 1
python3 ./PPI_fctPSN_overlap_calculations.py HI-III SEEK_matrix SeekPSN 100 PPI -1 2.1
python3 ./PPI_fctPSN_overlap_calculations.py HI-III seq_similarity_matrix SeqIdent 100 PPI 0 1
python3 ./PPI_fctPSN_overlap_calculations.py HI-III_PSN avana_2017_dep_corr_ENSG_ID coFitnessPSN 100 Jaccard_similarity 0 1 &
python3 ./PPI_fctPSN_overlap_calculations.py HI-III_PSN SEEK_matrix SeekPSN 100 Jaccard_similarity -1 2.1 &
python3 ./PPI_fctPSN_overlap_calculations.py HI-III_PSN seq_similarity_matrix SeqIdent 100 Jaccard_similarity 0 1 &
# produce all functional overlap and PSN related plots
jupyter nbconvert --to notebook --execute fct_overlap_analysis_figures.ipynb
rm fct_overlap_analysis_figures.nbconvert.ipynb


"""
code for analyses from integrating HuRI with structural data
"""
# code that draws (not computes) plots on distances of interaction interfaces to protein termini
jupyter nbconvert --to notebook --execute distance_terminus_interface_plots.ipynb
rm distance_terminus_interface_plots.nbconvert.ipynb


"""
code to draw heatmaps for network bias - correlation with gene property analysis and NCK2 brain expression heatmap
"""
jupyter nbconvert --to notebook --execute heatmaps_corr_NCK2.ipynb
rm heatmaps_corr_NCK2.nbconvert.ipynb


"""
code to generate tissue networks and perform analyses on it
"""
# generate tissue networks
python3 generate_TSNs.py
# get collapsed tissue networks for each source network -> all edges between proteins expressed together
# in at least one tissue
# including testis
python3 ./generate_CSN_pickles.py Lit-BM-17 True
python3 ./generate_CSN_pickles.py BioPlex True
python3 ./generate_CSN_pickles.py CoFrac True
python3 ./generate_CSN_pickles.py QUBIC True
python3 ./generate_CSN_pickles.py HI-III True
python3 ./generate_CSN_pickles.py HI-II-14 True
python3 ./generate_CSN_pickles.py HI-I-05 True
# excluding testis
python3 ./generate_CSN_pickles.py Lit-BM-17 False
python3 ./generate_CSN_pickles.py BioPlex False
python3 ./generate_CSN_pickles.py CoFrac False
python3 ./generate_CSN_pickles.py QUBIC False
python3 ./generate_CSN_pickles.py HI-III False
python3 ./generate_CSN_pickles.py HI-II-14 False
python3 ./generate_CSN_pickles.py HI-I-05 False
# generate random graphs
python3 ./random_graphs.py HI-III 1000 ../data/katjas_data/GTEx/analysis_GTExV6_2018_frozen/ collapsed_CSN_rand_networks 2 True
python3 ./random_graphs.py HI-III 1000 ../data/katjas_data/GTEx/analysis_GTExV6_2018_frozen/ collapsed_CSN_rand_networks_no_testis 2 False
python3 ./random_graphs.py HI-III 1000 ../data/katjas_data/GTEx/analysis_GTExV6_2018_frozen/ random_CSN_networks 1
python3 ./random_graphs.py HI-III 1000 ../data/katjas_data/GTEx/analysis_GTExV6_2018_frozen/ rand_source_networks 0
# calculate significances of TiP genes being close in tissue networks
python3 ./TiP_network_closeness.py HI-III nl,asp 2 False 1000
# produce all plots related to tissue network analyses
jupyter nbconvert --to notebook --execute tissue_network_analyses.ipynb
rm tissue_network_analyses.nbconvert.ipynb
# generate graph objects to draw TiP networks in cytoscape
python3 draw_TiP_networks.py HI-III all TiP_networks
# draw TiP networks with cytoscape
jupyter nbconvert --to notebook --execute draw_TiP_networks.ipynb
rm draw_TiP_networks.nbconvert.ipynb

"""
code to generate plots for analyses of apoptosis candidate genes and cell death assay data
"""
jupyter nbconvert --to notebook --execute apoptosis_candidate_plots.ipynb
rm apoptosis_candidate_plots.nbconvert.ipynb

"""
code to generate plots for analyses of tissue-specific Mendelian diseases
"""
# draw plots on stats of tissue-specific diseases and expression heatmaps
jupyter nbconvert --to notebook --execute tissue_spec_diseases.ipynb
rm tissue_spec_diseases.nbconvert.ipynb
# produce summary table with annotated edgotyping data
jupyter nbconvert --to notebook --execute get_HuRI_paper_edgotyping_data.ipynb
rm get_HuRI_paper_edgotyping_data.nbconvert.ipynb
# generate causal protein - TiP protein networks per tissue
jupyter nbconvert --to notebook --execute visualize_CG_TiP_data_by_tissue.ipynb
rm visualize_CG_TiP_data_by_tissue.nbconvert.ipynb
# compute significance of causal proteins to interact with TiP proteins
python3 ./test_disease_TiP_connectivity.py HI-III 1000 2
# plot results from closeness analysis of causal proteins to TiP proteins
jupyter nbconvert --to notebook --execute causal_gene_TiP_PPI_analysis.ipynb
rm causal_gene_TiP_PPI_analysis.nbconvert.ipynb
