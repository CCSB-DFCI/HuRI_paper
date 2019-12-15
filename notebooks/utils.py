import pandas as pd

# TODO
# huri.load_id_map('hgnc_id', 'ensembl_gene_id')
# huri.load_id_map('hgnc_symbol', 'ensembl_gene_id')


def load_nw_gsm():
    hi = pd.read_csv('../data/processed/Supplementary Table 4.txt', sep='\t')
    hi = hi.join(pd.get_dummies(hi['assay_id'],
                                prefix='in_assay_v', prefix_sep=''))
    for col in [c for c in hi.columns if c.startswith('in_assay_v')]:
        hi[col] = hi[col].astype(bool)
    for assay_id in hi['assay_id'].unique():
        nScreens = int(hi.loc[hi['assay_id'] == assay_id, 'screens']
                         .str.split(',', expand=True)
                         .astype(float).max().max())
        for i in map(str, range(1, nScreens + 1)):
            column = 'in_screen_v' + str(assay_id) + '_' + i
            hi[column] = (hi.loc[hi['assay_id'] == assay_id, 'screens']
                            .str.split(',')
                            .apply(lambda x: i in x))
            hi[column] = hi[column].fillna(False)
    hi = (hi.groupby(['ensembl_gene_id_a', 'ensembl_gene_id_b'])
            .agg({i: 'sum' for i in hi.columns if i.startswith('in_')})
            .reset_index())
    return hi


def load_nw_hi_union():
    df = pd.read_csv('../data/processed/Supplementary Table 9.txt', sep='\t')
    for column in df.columns:
        if column.startswith('in_'):
            df[column] = df[column] == 1
    df = df.rename(columns={'Ensembl_gene_id_a': 'ensembl_gene_id_a',
                            'Ensembl_gene_id_b': 'ensembl_gene_id_b'})
    return df


def load_nw_hi_iii():
    df = pd.read_csv('../data/processed/Supplementary Table 6.txt', sep='\t')
    for column in df.columns:
        if column.startswith('in_'):
            df[column] = df[column] == 1
    df = df.rename(columns={'Ensembl_gene_id_a': 'ensembl_gene_id_a',
                            'Ensembl_gene_id_b': 'ensembl_gene_id_b'})
    return df


def load_nw_lit_bm_17():
    return pd.read_csv('../data/processed/Supplementary Table 26.txt', sep='\t')


def load_space_iii(id_type='ensembl_gene_id'):
    if id_type == 'ensembl_gene_id':
        return set(pd.read_csv('../data/processed/Supplementary Table 20.txt', sep='\t')
                       ['ensembl_gene_id']
                       .str.replace(r'\..*', '')
                       .values)
    elif id_type == 'orf_id':
        return set(pd.read_csv('../data/processed/Supplementary Table 20.txt', sep='\t')
                     ['orf_id']
                     .values)
    else:
        raise ValueError('Unrecognized id_type')
        
        
def load_protein_coding_genome():       
    return set(pd.read_csv('../data/processed/Supplementary Table 1.txt', sep='\t')
                 ['Ensembl_gene_id'].values)


def load_subcellular_location(min_reliability='Approved'):
    """Subcellular location from Human Protein Atlas

    Args:
        min_reliability ({'Enhanced', 'Supported', 'Approved', 'Uncertain'}): minimum
            threshold to consider. For detailed description see:
            https://www.proteinatlas.org/about/assays+annotation

    Returns:
        DataFrame: indexed by ensembl_gene_id, boolean columns for whether
        protein found in compartment for every celular compartment.

    """
    reliabilities = ['Enhanced', 'Supported', 'Approved', 'Uncertain']
    if min_reliability not in reliabilities:
        raise ValueError(min_reliability + 'not valid value.\n' +
                         'Valid options are: ' + ','.join(reliabilities))
    file_path = '../data/external/subcellular_location.tsv'
    df = pd.read_table(file_path)
    to_use = reliabilities[:reliabilities.index(min_reliability) + 1]
    df['locations'] = df[to_use].apply(lambda x: ';'.join([i for i in x
                                                           if pd.notnull(i)]),
                                       axis=1)
    df = df.loc[df['locations'] != '', :]
    all_locations = set(df['locations'].str.split(';').sum())
    for loc in all_locations:
        df['in_' + loc] = df['locations'].apply(lambda x: loc in x.split(';'))
    df = (df.loc[:, ['Gene'] + ['in_' + i for i in all_locations]]
            .rename(columns={'Gene': 'ensembl_gene_id'}))
    if df['ensembl_gene_id'].duplicated().any():
        raise UserWarning('Unexpected duplicate ensembl gene IDs')
    df = df.set_index('ensembl_gene_id')
    return df


def load_number_publications_per_gene(protein_coding_only=True):
    """Number of unique PubMed IDs associated with each ensembl gene ID.

    Note:
        There can be some circularity when making samograms using this count
        since there are papers such as Rolland et al. here, as well as MGC
        papers and BioPlex, which uses the hORFeome. These can result in a
        small sparse zone for very unstudied genes.

    Args:
        protein_coding_only (bool): restrict to only protein coding genes

    Returns:
        pandas.Series: number of publications indexed by ensembl_gene_id

    """
    file_path = '../data/external/gene2pubmed'
    gn2pm = pd.read_table(file_path)
    gn2pm = gn2pm.loc[gn2pm['#tax_id'] == 9606, ['GeneID', 'PubMed_ID']]
    id_map = pd.read_csv('../data/external/gene2ensembl', sep='\t')
    id_map = (id_map.loc[id_map['#tax_id'] == 9606,
                         ['GeneID', 'Ensembl_gene_identifier']]
                    .drop_duplicates())
    gn2pm = pd.merge(gn2pm, id_map, how='inner', on='GeneID')
    gn2pm = gn2pm.drop(columns='GeneID').drop_duplicates()
    n_pubs = gn2pm.groupby('Ensembl_gene_identifier').size()
    pcg = load_protein_coding_genome()
    missing_vals = pcg.difference(set(n_pubs.index))
    n_pubs = pd.concat([n_pubs, pd.Series(index=list(missing_vals), data=0)])
    n_pubs.name = 'n_pubmed_ids'
    n_pubs.index.name = 'ensembl_gene_id'
    if protein_coding_only:
        n_pubs = n_pubs[n_pubs.index.isin(load_protein_coding_genome())]
    return n_pubs
