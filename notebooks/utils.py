from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import io
import os
import warnings
import hashlib
import base64
import glob

import numpy as np
from scipy import stats
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
import pandas as pd
from bioservices import BioMart
import igraph
from tqdm import tqdm


def generate_random_networks(real_nw, n=1000, cache=True):
    """Generate a set of degree-preserved randomized networks.
    The networks produced will be saved for re-use later.
    Warning:
        The produced networks will be a single fully connected component.
        You should check that these will be representative of your input
        network.
    Args:
        real_nw (igraph.Graph): Input network to randomize. Must have no
                                self-interactions
        n (int): Number of random networks to generate
        cache (bool): Controls whether to save the networks for re-use later
    Returns:
        list(igraph.Graph): randomly generated networks
    """
    if any(real_nw.is_loop()):
        msg = """real_nw must not contain self-interactions.
                 Remove with:
                 g.simplify()
                 g.delete_vertices([v for v in g.vs if g.degree() == 0])"""
        raise ValueError(msg)
    if any([d == 0 for d in real_nw.degree()]):
        msg = """real_nw must not contain any unconnected nodes.
                 Remove with:
                 g.delete_vertices([v for v in g.vs if v.degree() == 0])"""
        raise ValueError(msg)
    cache_dir = '../random_nw_cache'
    degree_seq = sorted(real_nw.degree())
    ds_hash = base64.urlsafe_b64encode(hashlib.md5(repr(degree_seq).encode('utf8')).digest()).decode('utf8')
    rand_nws = []
    i = 1
    while True:
        cache_dir = os.path.join(cache_dir,
                                 'random_networks',
                                 ds_hash + '_' + str(i))
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
            break
        else:
            rand_nws = [igraph.load(f) for f in glob.glob(cache_dir + '/random_network_*.txt')]
            if all([sorted(nw.degree()) == degree_seq for nw in rand_nws]):
                break
            else:
                i += 1
    if len(rand_nws) < n:
        for __ in tqdm(range(n - len(rand_nws))):
            rand_nws.append(igraph.Graph.Degree_Sequence(real_nw.degree(),
                                                         method='vl'))
    if cache:
        for i, rand_nw in enumerate(rand_nws):
            outpath = os.path.join(cache_dir,
                                   'random_network_' + str(i) + '.txt')
            rand_nw.write_edgelist(outpath)
    if all([real_nw.degree() == rand_nw.degree() for rand_nw in rand_nws]):
        for attribute in real_nw.vertex_attributes():
            for rand_nw in rand_nws:
                rand_nw.vs[attribute] = real_nw.vs[attribute]
    else:
        warnings.warn('WARNING: random network nodes in different order to real network')

        def get_degree(v):
            return v.degree()

        for attribute in real_nw.vertex_attributes():
            for rand_nw in rand_nws:
                for v_real, v_rand in zip(sorted(real_nw.vs, key=get_degree),
                                          sorted(rand_nw.vs, key=get_degree)):
                    v_rand[attribute] = v_real[attribute]
    return rand_nws

def merge_node_and_edge_tables(nodes,
                               edges,
                               id_type=None,
                               suffixes=('_a', '_b'),
                               node_id_column=None):
    """Combine data on nodes and edges into a table of edges.
    Args:
        nodes (pandas.DataFrame): table of nodes
        edges (pandas.DataFrame): table of edges
    """
    if id_type is None:
        _guess_node_id(edges.columns, suffixes)
    if node_id_column is None:
        df = pd.merge(edges,
                      nodes,
                      right_index=True,
                      left_on=id_type + suffixes[0],
                      how='left')
        df = pd.merge(df,
                      nodes,
                      right_index=True,
                      left_on=id_type + suffixes[1],
                      how='left',
                      suffixes=suffixes)
    else:
        df = pd.merge(edges,
                      nodes,
                      right_on=node_id_column,
                      left_on=id_type + suffixes[0],
                      how='left')
        df = pd.merge(df,
                      nodes,
                      right_on=node_id_column,
                      left_on=id_type + suffixes[1],
                      how='left',
                      suffixes=suffixes)
    if df.columns.duplicated().any():
        df = df.loc[:, ~df.columns.duplicated()]
    return df


def split_node_and_edge_tables(edges, id_type=None, suffixes=None):
    es = edges.copy()
    if id_type is None:
        if suffixes is None:
            id_type, suffixes = _guess_id_type_and_suffixes(es.columns)
        else:
            id_type = _guess_node_id(es.columns, suffixes)
    elif suffixes is None:
        suffixes = es.columns[0].replace(id_type, ''), es.columns[1].replace(id_type, '')
    columns_a = [c for c in es.columns if c.endswith(suffixes[0])]
    columns_b = [c for c in es.columns if c.endswith(suffixes[1])]    
    ns_a = (es.loc[:, columns_a]
              .rename(columns={c: c[:-len(suffixes[0])] for c in columns_a})
              .drop_duplicates())
    ns_b = (es.loc[:, columns_b]
              .rename(columns={c: c[:-len(suffixes[1])] for c in columns_b})
              .drop_duplicates())
    ns = (pd.concat([ns_a, ns_b])
            .drop_duplicates()
            .set_index(id_type, verify_integrity=True))
    gene_id_columns = [id_type + s for s in suffixes]
    es = es.loc[:, gene_id_columns + [c for c in es.columns if c not in (columns_a + columns_b)]]
    return ns, es


def _guess_node_id(columns, suffixes=('_a', '_b')):
    """
    Given a table of edges, try and guess the node IDs.
    Should I give a warning? i.e. guessing gene id type is x, to silence pass id_type=x?
    Maybe only give warning if it's ambiguous? i.e. there are not multiple?
    """
    cols_a = [c[:-len(suffixes[0])] for c in columns if c.endswith(suffixes[0])]
    cols_b = [c[:-len(suffixes[1])] for c in columns if c.endswith(suffixes[1])]
    for col_a in cols_a:
        for col_b in cols_b:
            if col_a == col_b:
                return col_a
    raise UserWarning('Could not guess node IDs from: ' + ' | '.join(columns))


def _guess_id_type_and_suffixes(columns):
    id_a, id_b = columns[:2]
    if id_a.split('_')[:-1] == id_b.split('_')[:-1]:
        return ('_'.join(id_a.split('_')[:-1]),
                ('_' + id_a.split('_')[-1],
                 '_' + id_b.split('_')[-1]))
    else:
        raise UserWarning('Failed to guess id_type and suffixes. Please specify.')


def format_network(data, fmt, id_type=None, suffixes=('_a', '_b'), directed=False):
    """Convert data format of network between pandas/networkx/igraph etc.
    fmt options are:
        - `pandas`: DataFrame
        - `igraph`: igraph Graph
        - `list`: list of a dict for each PPI
    Args:
        data (pandas.DataFrame / igraph.Graph / list): input network
        fmt (str): format to convert to; pandas/nx/igraph/list
        id_type (str): gene/protein identifier type used
        suffixes (tuple(str, str)): at the end of id_type to distinguish the two nodes
        directed (bool): return directed for igraph/networkx
    Returns:
        (pandas.DataFrame / igraph.Graph / list): network in specified format
    """
    valid_fmt = ['pandas', 'igraph', 'list']
    if fmt not in valid_fmt:
        raise UserWarning('Unsupported fmt: ' + fmt +
                          '\nValid options are: ' + '/'.join(valid_fmt))
    fmt_out = fmt
    fmt_in = 'unknown'
    for fmt_desc, data_type in [('pandas', pd.DataFrame),
                                ('igraph', igraph.Graph),
                                ('list', list)]:
        if isinstance(data, data_type):
            fmt_in = fmt_desc
            break
    if fmt_in == 'unknown':
        raise ValueError('Unsupported input type: ' + type(data))
    if id_type is None:
        if fmt_in == 'pandas':
            id_type, suffixes = _guess_id_type_and_suffixes(data.columns)
        else:
            raise ValueError('Require value for id_type argument')
    id_a = id_type + suffixes[0]
    id_b = id_type + suffixes[1]

    if fmt_in != 'pandas' and fmt_out != 'pandas':
        # via pandas.DataFrame, so don't have to code every possible conversion
        tbl = format_network(data,
                             'pandas',
                             id_type=id_type,
                             suffixes=suffixes,
                             directed=directed)
        return format_network(tbl,
                              fmt_out,
                              id_type=id_type,
                              suffixes=suffixes,
                              directed=directed)

    elif fmt_in == 'pandas' and fmt_out == 'pandas':
        return data

    elif fmt_in == 'pandas' and fmt_out == 'igraph':
        node_df, edge_df = split_node_and_edge_tables(data, id_type=id_type, suffixes=suffixes)
        g = igraph.Graph()
        g = g.TupleList([(a, b) for a, b in edge_df[[id_a, id_b]].values],
                        directed=directed)
        for column in edge_df.columns:
            if column not in [id_a, id_b]:
                g.es[column] = edge_df[column].values
        for column in node_df.columns:
            g.vs[column] = node_df.loc[g.vs['name'], column].values
        return g

    elif fmt_in == 'pandas' and fmt_out == 'list':
        d = data.to_dict()
        mylist = [{k: d[k][idx] for k in d.keys()} for idx in d[id_a].keys()]
        return mylist

    elif fmt_in == 'igraph' and fmt_out == 'pandas':
        d = {id_a: [data.vs[e.source]['name'] for e in data.es],
             id_b: [data.vs[e.target]['name'] for e in data.es]}
        d.update({k: data.es[k] for k in data.es.attribute_names()})
        edge_df = pd.DataFrame(d)
        edge_df = edge_df.set_index(edge_df[id_a] + '_' + edge_df[id_b])
        node_df = pd.DataFrame({k: data.vs[k] for k in data.vs.attribute_names()}).set_index('name')
        out = merge_node_and_edge_tables(node_df, edge_df, id_type=id_type, suffixes=suffixes)
        return out

    elif fmt_in == 'list' and fmt_out == 'pandas':
        out = pd.DataFrame(data)
        out = out.set_index(out[id_a] + '_' + out[id_b])
        return out

    else:
        raise UserWarning('Something went wrong...')


def load_nw_qubic(id_type='ensembl_gene_id', fmt='pandas'):
    """QUBIC
    From [1]_
    Note:
        There are missing values in the stoichemetry columns in their data.
    References:
        .. [1] Hein et al., (2015). A Human Interactome in Three Quantitative
           Dimensions Organized by Stoichiometries and Abundances. Cell,
           163(3), 712-723. http://doi.org/10.1016/j.cell.2015.09.053
    Args:
        id_type ({'orf_id', 'ensembl_gene_id', 'uniprot_ac'}): protein/gene identifiers
        fmt ({'pandas', 'igraph', 'list'}): format to return network in
        keep_bait_prey (bool): retain the bait prey information
    Returns:
        :class:`pandas.DataFrame` / igraph.Graph / list:
        QUBIC - one entry per pair

    """
    if id_type not in ['ensembl_gene_id']:
        raise ValueError('Unsupported ID type: ' + id_type)
    suffixes = ('_a', '_b')
    fpath = ('../data/external/Hein_et_al_Cell_2015_TableS2.xlsx')
    qb = pd.read_excel(fpath,
                       sheet_name='interactions',
                       usecols=[2, 3, 11, 12])
    qb = (qb.join(qb['bait.IDs'].str.split(';', expand=True)
                                .stack()
                                .rename('uniprot_ac_bait')
                                .reset_index(level=1, drop=True))
            .drop(columns='bait.IDs'))
    qb = (qb.join(qb['prey.IDs'].str.split(';', expand=True)
                                .stack()
                                .rename('uniprot_ac_prey')
                                .reset_index(level=1, drop=True))
            .drop(columns='prey.IDs'))
    qb = qb.reset_index(drop=True)
    qb['uniprot_ac_bait'] = qb['uniprot_ac_bait'].str.replace('-.*', '')
    qb['uniprot_ac_prey'] = qb['uniprot_ac_prey'].str.replace('-.*', '')
    qb['uniprot_ac_a'] = (qb.loc[:, ['uniprot_ac_bait',
                                     'uniprot_ac_prey']]
                            .min(axis=1))
    qb['uniprot_ac_b'] = (qb.loc[:, ['uniprot_ac_bait',
                                        'uniprot_ac_prey']]
                            .max(axis=1))
    qb = qb.drop(['uniprot_ac_bait', 'uniprot_ac_prey'], axis=1)
    ids = ['uniprot_ac_a', 'uniprot_ac_b']
    qb = qb.loc[:, ids + [c for c in qb.columns if c not in ids]]
    qb = qb.drop_duplicates(['uniprot_ac_a', 'uniprot_ac_b'])
    qb = qb.set_index(qb['uniprot_ac_a'] + '_' + qb['uniprot_ac_b'])
    qb = qb.drop_duplicates()
    if id_type == 'ensembl_gene_id':
        qb = map_nw_ids(qb, 'uniprot_ac', 'ensembl_gene_id',
                        suffixes=suffixes, directed=False)
        pcg = load_protein_coding_genome()
        qb = qb.loc[qb['ensembl_gene_id' + suffixes[0]].isin(pcg) &
                    qb['ensembl_gene_id' + suffixes[1]].isin(pcg), :]
    qb = format_network(qb,
                        fmt,
                        id_type=id_type,
                        suffixes=suffixes,
                        directed=False)
    return qb


def load_nw_cofrac(id_type='ensembl_gene_id', fmt='pandas'):
    """Interactome based on co-fractionation and machine learning.

    References:
        .. [1] Wan et al., (2015). Panorama of ancient metazoan macromolecular
           complexes. Nature, 525(7569), 339-344.
           http://doi.org/10.1038/nature14877
    Args:
        id_type ({'orf_id', 'ensembl_gene_id'}): protein/gene identifiers
        fmt ({'pandas', 'nx', 'igraph', 'list'}): format to return network in
    Returns:
        :class:`pandas.DataFrame` / igraph.Graph / list:
        CoFrac - one entry per pair
    """
    if id_type not in ['ensembl_gene_id']:
        raise ValueError('Unsupported ID type: ' + id_type)
    fpath = '../data/external/High_confidence_16655_correlations_and_ppi_scores.txt'
    cf = (pd.read_csv(fpath, usecols=['id1', 'id2'], sep='\t')
            .rename(columns={'id1': 'ensembl_gene_id_a',
                                'id2': 'ensembl_gene_id_b'}))
    a = cf.loc[:, ['ensembl_gene_id_a', 'ensembl_gene_id_b']].min(axis=1)
    b = cf.loc[:, ['ensembl_gene_id_a', 'ensembl_gene_id_b']].max(axis=1)
    cf['ensembl_gene_id_a'] = a
    cf['ensembl_gene_id_b'] = b
    cf = cf.drop_duplicates()
    cf = cf.set_index(cf['ensembl_gene_id_a'] +
                        '_' +
                        cf['ensembl_gene_id_b'])
    if id_type == 'orf_id':
        cf = map_nw_ids(cf, 'ensembl_gene_id', 'orf_id')
    elif id_type == 'uniprot_ac':
        cf = map_nw_ids(cf, 'ensembl_gene_id', 'uniprot_id')
    elif id_type == 'ensembl_gene_id':
        pcg = load_protein_coding_genome()
        cf = cf.loc[cf['ensembl_gene_id_a'].isin(pcg) &
                    cf['ensembl_gene_id_b'].isin(pcg), :]
    cf = format_network(cf, fmt, id_type=id_type)
    return cf


def load_nw_bioplex(id_type='ensembl_gene_id', fmt='pandas'):
    """The BioPlex 2.0 AP-MS interactome.
    This dataset is from [1]_ which includes a reanlysis of the original
    paper [2]_ and does not include the unpublished data.
    References:
        .. [1] Huttlin et al., (2017). Architecture of the human interactome
           defines protein communities and disease networks. Nature,
           545(7655), 505-509. http://doi.org/10.1038/nature22366
        .. [2] Huttlin et al., (2015). The BioPlex Network: A Systematic
           Exploration of the Human Interactome. Cell, 162(2), 425-440.
           http://doi.org/10.1016/j.cell.2015.06.043
    Note:
        The original BioPlex file contains no homodimers, however, after
        mapping to ensembl IDs / ORF IDs there are a very small number of
        self-interactions due to inconsistencies between different gene
        annotations in a highly variable, hard-to-annotate genomic region.
    Args:
        id_type ({'orf_id', 'ensembl_gene_id', 'uniprot_ac'}): protein/gene identifiers
        fmt ({'pandas', 'nx', 'igraph', 'list'}): format to return network in
        version ({'1.0', '2.0', 'unpublished', '2.0_unfiltered'}): available
            versions are the 2015 and 2017 published versions, the unpublished
            data and the unfiltered version of 2.0.
        keep_bait_prey (bool): retain the bait prey information
    Returns:
        :class:`pandas.DataFrame` / igraph.Graph /:class:`networkx.Graph` / list:
        BioPlex 2.0 - one entry per pair

    """
    if id_type not in ['ensembl_gene_id']:
        raise ValueError('Unsupported ID type: ' + id_type)
    fpath = '../data/external/BioPlex_interactionList_v4a.tsv'
    suffixes = ('_a', '_b')

    # Learnt that A is bait and B is prey in this file from email
    # communication between Katja and Ed Huttlin. Caveat is that
    # the pairs are unique in that file so in cases where both ways
    # are found (bait-X/prey-Y & bait-Y/prey-X), then only one is chosen.
    col_bait = 'UniprotA'
    col_prey = 'UniprotB'
    bp = (pd.read_csv(fpath,
                      sep='\t',
                      header=0,
                      usecols=[col_bait, col_prey])
            .rename(columns={col_bait: 'uniprot_ac_bait',
                             col_prey: 'uniprot_ac_prey'}))
    bp = bp.loc[(bp['uniprot_ac_bait'] != 'UNKNOWN') &
                (bp['uniprot_ac_prey'] != 'UNKNOWN'), :]
    # converting from isoform level UniProt AC to protein/gene level
    bp['uniprot_ac_bait'] = bp['uniprot_ac_bait'].str.replace('-.*', '')
    bp['uniprot_ac_prey'] = bp['uniprot_ac_prey'].str.replace('-.*', '')
    bp['uniprot_ac_a'] = (bp.loc[:, ['uniprot_ac_bait',
                                        'uniprot_ac_prey']]
                            .min(axis=1))
    bp['uniprot_ac_b'] = (bp.loc[:, ['uniprot_ac_bait',
                                        'uniprot_ac_prey']]
                            .max(axis=1))
    bp = bp.drop(['uniprot_ac_bait', 'uniprot_ac_prey'], axis=1)
    bp = bp.drop_duplicates()
    bp = bp.set_index(bp['uniprot_ac_a'] + '_' + bp['uniprot_ac_b'])
    bp = bp.drop_duplicates()
    if id_type == 'orf_id':
        bp = map_nw_ids(bp, 'uniprot_ac', 'orf_id', via='ensembl_gene_id',
                        suffixes=suffixes, directed=False)
    elif id_type == 'ensembl_gene_id':
        bp = map_nw_ids(bp, 'uniprot_ac', 'ensembl_gene_id',
                        suffixes=suffixes, directed=False)
        pcg = load_protein_coding_genome()
        bp = bp.loc[bp['ensembl_gene_id' + suffixes[0]].isin(pcg) &
                    bp['ensembl_gene_id' + suffixes[1]].isin(pcg), :]

    bp = format_network(bp,
                        fmt,
                        id_type=id_type,
                        suffixes=suffixes,
                        directed=False)
    return bp


def map_nw_ids(df, id_in, id_out, via=None, suffixes=('_a', '_b'),
               directed=False, space_iii=True, agg=None):
    """Tables mapping between different protein/ORF/gene identifiers.
    All combinations of pairs of IDs are returned. Mappings are symmetric.
    Supported IDs are: orf_id, ensembl_gene_id, uniprot_ac
    Where orf_id refers to our internal ORF IDs.
    The ORF ID mappings come from the ORFeome annotation project where ORFs
    were sequence aligned to ensembl transcripts and proteins. The uniprot
    to ensembl mappings are from a file provided by UniProt.
    Note:
        The ensembl gene IDs and UniProt ACs are an unfiltered redundant set.
        You will probably want to filter the set after mapping.
    Args:
        df (DataFrame): table of protein-protein interactions,
                        one row for each pair
        id_in/id_out (str): gene/protein identifiers to map between
        via (str): optional identifier to use as an intermediate in the
                   mapping. Useful for example when mapping uniprot IDs to
                   our ORF IDs from experiments that can ony determine the
                   gene-level then can map via ensembl_gene_id.
        suffixes (tuple(str, str)): of the two protein identifiers of each pair
        directed (bool): if True, retain the mapping of the two proteins,
                         if False, sort the IDs and drop duplicates
        space_iii (bool): restrict ORF IDs to Space III, where there is one
                             ORF per gene. Only has an effect if one of id_in,
                             id_out or via_ID is 'orf_id'.
        agg (function): Optional. Custom aggregation function to choose a
                        single pair or combine multiple pairs when input pairs
                        are mapped to the same pair in the new ID. Function
                        must take a DataFrame and return a DataFrame with no
                        duplicate pairs.
    Returns:
        DataFrame: PPI dataset mapped to new protein/gene ID. All unique
                   combinations of the new ID are mapped to, so one pair in
                   the input dataset can map to multiple pairs in the output
                   and vice-versa.
    """
    id_in_a = id_in + suffixes[0]
    id_in_b = id_in + suffixes[1]
    id_out_a = id_out + suffixes[0]
    id_out_b = id_out + suffixes[1]
    id_map = load_id_map(id_in, id_out)
    out = df.copy()
    out = (pd.merge(out, id_map,
                    how='inner',
                    left_on=id_in_a,
                    right_on=id_in)
             .drop(columns=id_in)
             .rename(columns={id_out: id_out_a}))
    out = (pd.merge(out, id_map,
                    how='inner',
                    left_on=id_in_b,
                    right_on=id_in)
             .drop(columns=id_in)
             .rename(columns={id_out: id_out_b}))
    if out.loc[:, [id_out_a, id_out_b]].isnull().any().any():
        raise UserWarning('Unexpected missing values')
    if not directed:
        a = out[[id_out_a, id_out_b]].min(axis=1)
        b = out[[id_out_a, id_out_b]].max(axis=1)
        out[id_out_a] = a
        out[id_out_b] = b
    out = out.drop(columns=[id_in_a, id_in_b])
    pair_duplicates = out.duplicated(subset=[id_out_a, id_out_b], keep=False)
    if (out.duplicated(keep=False) != pair_duplicates).any() and agg is None:
        warnings.warn('Warning: mapping between gene/protein identifiers has '
                      'resulted in different pairs in the input ID being mapped to '
                      'the same pair in the output ID.\n'
                      'You may wish to use the `agg` argument to customize '
                      'the choice of which of the pair\'s infomation to keep or how '
                      'to combine the information from multiple pairs.')
    if agg is None:
        out = out.drop_duplicates(subset=[id_out_a, id_out_b])
    else:
        out = agg(out)
        if out.duplicated().any():
            raise ValueError('Problem with your agg function')
    cols = out.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    out = out.loc[:, cols]
    out = out.set_index(out[id_out_a].astype(str) +
                        '_' +
                        out[id_out_b].astype(str))
    return out


def ensembl_human_id_mappings(release_date='aug2017'):
    bm = BioMart(host=release_date + '.archive.ensembl.org')
    bm.datasets('ENSEMBL_MART_ENSEMBL')
    bm.add_dataset_to_xml('hsapiens_gene_ensembl')
    bm.add_attribute_to_xml('ensembl_gene_id')
    bm.add_attribute_to_xml('ensembl_transcript_id')
    # can have maximum 3 external identifiers
    bm.add_attribute_to_xml('hgnc_id')
    bm.add_attribute_to_xml('hgnc_symbol')
    if 'entrezgene' in bm.attributes(dataset='hsapiens_gene_ensembl').keys():
        bm.add_attribute_to_xml('entrezgene')
    elif 'entrezgene_id' in bm.attributes(dataset='hsapiens_gene_ensembl').keys():
        bm.add_attribute_to_xml('entrezgene_id')
    else:
        raise UserWarning('couldnt find correct name for entrez gene id field')
    res = bm.query(bm.get_xml())
    df = pd.read_csv(io.StringIO(res),
                     names=['ensembl_gene_id',
                            'ensembl_transcript_id',
                            'hgnc_id',
                            'hgnc_symbol',
                            'entrez_gene_id'],
                     sep='\t')
    return df


def load_id_map(id_in, id_out):
    valid_ids = ['ensembl_gene_id',
                 'uniprot_ac',
                 'hgnc_id',
                 'hgnc_symbol']
    for id_type in [id_in, id_out]:
        if id_type not in valid_ids:
            raise ValueError('Unsupported ID: ' + id_type +
                             '\nChoices are: ' + ', '.join(valid_ids))
    if id_in == id_out:
        raise ValueError('Invalid arguments: id_in == id_out')
    ids = set([id_in, id_out])
    if ids == {'uniprot_ac', 'ensembl_gene_id'}:
        return pd.read_csv('../data/external/uniprot_to_ensembl.tsv', sep='\t')
    id_map_path = '../data/external/ensembl_v90_gene_id_mapping.tsv'
    if os.path.exists(id_map_path):
        df = pd.read_csv(id_map_path, sep='\t')
    else:
        df = ensembl_human_id_mappings()
        df.to_csv(id_map_path, sep='\t', index=False)
    df = df.loc[:, [id_in, id_out]].dropna().drop_duplicates()
    if id_in == 'hgnc_id' or id_out == 'hgnc_id':
        df['hgnc_id'] = df['hgnc_id'].str.replace('HGNC:', '').astype(int)
    if id_in == 'entrez_gene_id' or id_out == 'entrez_gene_id':
        df['entrez_gene_id'] = df['entrez_gene_id'].astype(int)
    return df


def load_nw_gsm():
    hi = pd.read_csv('../data/processed/Supplementary Table 7.txt', sep='\t')
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


def load_nw_hi_union(id_type='ensembl_gene_id', fmt='pandas'):
    df = pd.read_csv('../data/processed/Supplementary Table 11.txt', sep='\t')
    for column in df.columns:
        if column.startswith('in_'):
            df[column] = df[column] == 1
    df = df.rename(columns={'Ensembl_gene_id_a': 'ensembl_gene_id_a',
                            'Ensembl_gene_id_b': 'ensembl_gene_id_b'})
    df = format_network(df,
                        fmt,
                        id_type=id_type,
                        suffixes=('_a', '_b'),
                        directed=False)
    return df


def load_nw_hi_iii(id_type='ensembl_gene_id', fmt='pandas'):
    df = pd.read_csv('../data/processed/Supplementary Table 9.txt', sep='\t')
    for column in df.columns:
        if column.startswith('in_'):
            df[column] = df[column] == 1
    df = df.rename(columns={'Ensembl_gene_id_a': 'ensembl_gene_id_a',
                            'Ensembl_gene_id_b': 'ensembl_gene_id_b'})
    df = format_network(df,
                        fmt,
                        id_type=id_type,
                        suffixes=('_a', '_b'),
                        directed=False)
    return df


def load_nw_lit_bm_17(id_type='ensembl_gene_id', fmt='pandas'):
    df = pd.read_csv('../data/processed/Supplementary Table 14.txt', sep='\t')
    df = format_network(df,
                        fmt,
                        id_type=id_type,
                        suffixes=('_a', '_b'),
                        directed=False)
    return df


def load_space_iii(id_type='ensembl_gene_id'):
    if id_type == 'ensembl_gene_id':
        return set(pd.read_csv('../data/processed/Supplementary Table 2.txt', sep='\t')
                       ['ensembl_gene_id']
                       .str.replace(r'\..*', '')
                       .values)
    elif id_type == 'orf_id':
        return set(pd.read_csv('../data/processed/Supplementary Table 2.txt', sep='\t')
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


def load_number_publications_per_gene():
    """Number of unique PubMed IDs associated with each ensembl gene ID.

    Note:
        There can be some circularity when making samograms using this count
        since there are papers such as Rolland et al. here, as well as MGC
        papers and BioPlex, which uses the hORFeome. These can result in a
        small sparse zone for very unstudied genes.

    Returns:
        pandas.Series: number of publications indexed by ensembl_gene_id

    """
    """
    # code used to make file from NCBI data:
    protein_coding_only = True
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
    """
    n_pubs = pd.read_csv('../data/processed/number_publications_per_gene.tsv',
                         index_col='ensembl_gene_id',
                         sep='\t',
                         squeeze=True)
    return n_pubs


def degree_per_protein(df, id_a='orf_id_a', id_b='orf_id_b'):
    """Calculate degree of each ORF
    Counts once for a self-interaction.
    Args:
        df (DataFrame): PPI data, one row per pair.
    Returns:
        Series: degree per protein
    """
    orfs = np.unique(np.concatenate([df[id_a].unique(),
                                     df[id_b].unique()]))
    d = pd.Series(index=orfs, data=np.zeros(orfs.shape))
    a = df.groupby(id_a).size()
    b = df.groupby(id_b).size()
    homo = df.loc[(df[id_a] == df[id_b]), id_a].values
    d.loc[a.index] = a
    d.loc[b.index] = d.loc[b.index] + b
    d.loc[homo] = d.loc[homo] - 1
    d.name = 'degree'
    return d


def validation_plot(positives=None,
                    n_tested=None,
                    data=None,
                    selections=None,
                    result_column='result',
                    labels=None,
                    colors=None,
                    ax=None,
                    bayes_errors=True,
                    y_max=1.0,
                    draw_numbers=True,
                    xlabel_rotation=0,
                    errorbar_capsize=0.9,
                    errorbar_thickness=None,
                    bar_spacing=0.2):
    """Compare the validation rate of different cateogies.
    
    Missing values are not used in the denominator.
    
    Args:
        positives (list(int)): number tested positive in each category
        n_tested (list(int)): number successfully tested in each category
        data (pandas.DataFrame): results of validation experiment
        selections (list(pandas.Series)): boolean rows index for each category
        ax (matplotlib.axes.Axes): Axes to draw plot onto
        bayes_errors (bool): do Bayesian error bars, if false use standard error
                            on proportion
        result_column (str): column containing 0/1/nan for result of test
        labels (list(str)): name of each category
        colors (list(str)): color for each bar
        y_max (float): y axis upper limit
        draw_numbers (bool): flag to print the numbers on top of the bars
        xlabel_rotation (float): rotating the x axis labels
        errorbar_capsize (float): as fraction of the width of the bars
        errorbar_thickness (float): width of error bar lines
        bar_spacing (float): must be between 0 and 1
    Examples:
        There are two ways to call the function. Either give it the raw data:
        .. plot::
            :context: close-figs
            >>> import ccsblib.ccsbplotlib as cplt
            >>> cplt.validation_plot(positives=[20, 1, 19],
            ...                      n_tested=[100, 100, 100],
            ...                      labels=['PRS', 'RRS', 'Y2H'],
            ...                      y_max=0.3)
        Or pass it the validation results as a DataFrame and a list of the rows
        for each category:
        .. plot::
            :context: close-figs
            >>> from ccsblib import huri
            >>> data = huri.load_validation_data()
            >>> exp = (data['assay'] == 'MAPPIT') & (data['standard_batch'] == 'Hvs01')
            >>> sources = ['lit_bm_2013_rand250', 'Hs01', 'RRS']
            >>> categories = [exp & (data['source'] == cat) for cat in sources]
            >>> cplt.validation_plot(data=data,
            ...                      selections=categories,
            ...                      labels=sources,
            ...                      y_max=0.25)
    """
    signature_a = positives is not None and n_tested is not None
    signature_b = data is not None and selections is not None
    if signature_a == signature_b:
        msg = """Must supply only one of both positives and n_tested or
                 both data and selections."""
        raise ValueError(msg)
    if signature_b:
        if not data.loc[data[result_column].notnull(), result_column].isin({0, 1}).all():
            raise ValueError('Only expect 0/1/missing in result column')
        positives = [(data.loc[rows, :][result_column] == 1).sum()
                     for rows in selections]
        n_tested = [data.loc[rows, :][result_column].notnull().sum()
                    for rows in selections]
    _validation_plot(positives=positives,
                     n_tested=n_tested,
                     colors=colors,
                     labels=labels,
                     ax=ax,
                     bayes_errors=bayes_errors,
                     y_max=y_max,
                     draw_numbers=draw_numbers,
                     xlabel_rotation=xlabel_rotation,
                     errorbar_capsize=errorbar_capsize,
                     errorbar_thickness=errorbar_thickness,
                     bar_spacing=bar_spacing)


def _validation_plot(positives,
                     n_tested,
                     colors=None,
                     labels=None,
                     ax=None,
                     bayes_errors=True,
                     y_max=1.0,
                     draw_numbers=True,
                     xlabel_rotation=0.,
                     errorbar_capsize=5,
                     errorbar_thickness=None,
                     bar_spacing=0.2):
    if len(positives) != len(n_tested):
        raise ValueError('Lengths of positives and n_tested must be equal')
    if any([p > n for p, n in zip(positives, n_tested)]):
        raise ValueError('Number of positives must be <= number tested')
    if bar_spacing > 1. or bar_spacing < 0.:
        msg = 'bar_spacing={}\nbar_spacing must be between 0 and 1'
        msg = msg.format(bar_spacing)
        raise ValueError(msg)
    bar_width = 1. - bar_spacing
    if ax is None:
        ax = plt.gca()
    if labels is None:
        labels = [''] * len(positives)
    if colors is None:
        colors = [None] * len(positives)
    ax.set_yticks(np.arange(0.0, 1.0, 0.1), minor=False)
    ax.set_yticks(np.arange(0.05, 1.0, 0.1), minor=True)
    ax.set_facecolor('0.96')
    ax.set_axisbelow(True)
    ax.grid(color='white', axis='y', which='both', zorder=5)
    pos = np.array(positives)
    tested = np.array(n_tested)
    neg = tested - pos
    fracs = pos / tested
    if bayes_errors:
        intv = stats.beta.interval(0.6827, pos + 1, neg + 1)
        errs = [fracs - intv[0], intv[1] - fracs]
        errs[0][pos == 0] = 0.
        errs[1][neg == 0] = 0.
    else:
        stdErrProp = np.sqrt((fracs * (1. - fracs)) / (pos + neg))
        errs = [stdErrProp, stdErrProp]
    for i in range(len(positives)):
        ax.bar(i,
               fracs[i],
               color=colors[i],
               label=labels[i],
               width=bar_width)
        if draw_numbers:
            c = 'white'
            h = 0.02  # default height to draw numbers
            if fracs[i] < h:
                c = 'black'
            if (errs[1][i] + fracs[i]) > h and (fracs[i] - errs[0][i]) < (h + 0.04):
                c = 'black'
                h = fracs[i] + errs[1][i] + 0.02
            ax.text(i, h, '{}/{}'.format(pos[i], pos[i] + neg[i]),
                    color=c, ha='center')
    bar_width_pixels = (ax.transData.transform((bar_width, 0)) -
                        ax.transData.transform((0, 0)))[0]
    # 72 comes from definition of a point as 1 / 72 inches
    bar_width_points = (72. / ax.figure.dpi) * bar_width_pixels
    ax.errorbar(range(fracs.shape[0]),
                fracs,
                yerr=errs,
                color='black',
                fmt='none',
                # I don't understand why I needed the factor of 0.5 below
                capsize=bar_width_points * errorbar_capsize * 0.5,
                elinewidth=errorbar_thickness,
                capthick=errorbar_thickness)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=xlabel_rotation)
    ax.set_ylim((0., y_max))
    ax.set_ylabel('Fraction positive')


def samogram(ppis, gene_property, n_bins,
             ax=None,
             id_a=None, id_b=None,
             correct_diagonal=False,
             self_interactions=True,
             draw_up=True, draw_right=False, draw_scale=True,
             reverse=False,
             vmax=None,
             ylim=None, ylabel=None, yticks=None,
             log=False, logy=False,
             zticks=None,
             cmap='Purples', color='purple',
             size_ratio=0.1, pad=0.04,
             colorbar_width=0.25,
             label=None, label_color=None):
    """2D histogram of PPIs with genes ordered by a given property
    Bar charts of the median of the property per bin can optionally be drawn.
    Gene property values are shuffled before sorting, so that genes with equal
    values are in a random order.
    See the :ref:`tutorial </samogram.ipynb>` for more information.
    Note:
        The diagonal bins have approximately half the number of possible
        interactions as the other bins, since the interactions are undirected.
        Specifically they have :math:`n^2/2 + n/2` possible interactions, with
        the other bins having :math:`n^2`, where :math:`n` is the number of
        genes per bin. This effect is counteracted by the much higher
        likelihood for genes to self-interact than interact with a
        random other gene. The smaller the bin size, the more
        self-interactions will have a visible impact.
    Args:
        ppis (pandas.DataFrame): table of interactions, one row per pair
        gene_property (pandas.Series): per-gene quantity to order proteins by
        n_bins (int): number of bins, heatmap will be n_bins x n_bins
        ax (matplotlib.axes.Axes): Axes to draw plot onto
        id_a/id_b (str): name of column containing gene/protein identifiers
        correct_diagonal (bool): increase number of interactions in diagonal
                                 bins to account for there being less possible
                                 combinations. Correction is dependent on value
                                 of `self_interactions` argument.
        self_interactions (bool): whether self interactions are plotted.
        draw_up/draw_right (bool): whether to add bar charts above / to the right
        draw_scale (bool): whether to draw the color scale
        reverse (bool): if true flips direction of ranking of genes
        vmax (int): upper limit on heatmap
        ylim ((float, float)): bar chart axis limits
        ylabel (str): bar chart axis label
        yticks/zticks (list(float)): bar chart / color scale axis ticks
        log/logy (bool): log scale for heatmap / bar charts
        cmap (str): colormap for heatmap
        color (str): color for bar charts
        size_ratio (float): fraction of size of bar charts to heatmap
        pad (float): space between bar charts and heatmap
        colorbar_width (float): size of colorbar as fraction of axes width
        label (str): label to put next to heatmap
    Returns:
        (list(matplotlib.axes.Axes), numpy.ndarray): new axes created, number of PPIs in each bin
    See Also:
        :func:`samogram_double`: plot two samograms in a single square
    Examples:
        Plot samogram of HI-III with default options:
        .. plot::
            :context: close-figs
            >>> import ccsblib.ccsbplotlib as cplt
            >>> from ccsblib import huri
            >>> hi_iii = huri.load_nw_hi_iii(id_type='ensembl_gene_id')
            >>> n_pub = huri.load_number_publications_per_gene()
            >>> cplt.samogram(hi_iii,
            ...               n_pub,
            ...               n_bins=40,
            ...               id_a='ensembl_gene_id_a',
            ...               id_b='ensembl_gene_id_b')
            >>> import matplotlib.pyplot as plt; fig = plt.gcf(); fig.set_size_inches(6, 6)
    """
    if size_ratio <= 0. or size_ratio >= 1.:
        raise ValueError('size_ratio must be between 0 and 1')
    if ax is None:
        ax = plt.gca()
    if id_a is None:
        id_a = gene_property.index.name + '_a'
    if id_b is None:
        id_b = gene_property.index.name + '_b'
    cmap = plt.get_cmap(cmap)
    cmap.set_under(color=(0, 0, 0, 0))  # fully transparent
    genes_ranked = (gene_property.sample(frac=1)  # randomize order first
                                 .sort_values()
                                 .to_frame())
    genes_ranked['rank'] = genes_ranked.reset_index().index.values
    binning = pd.qcut(genes_ranked['rank'], n_bins, labels=False)
    if reverse:
        binning = (n_bins - 1) - binning
    binned_ppis = (ppis.loc[ppis[id_a].isin(gene_property.index) &
                            ppis[id_b].isin(gene_property.index) &
                            ((ppis[id_a] != ppis[id_b]) | self_interactions),
                            [id_a, id_b]]
                       .apply(lambda x: x.map(binning), axis=0))
    binned_ppis_x = binned_ppis.min(axis=1)
    binned_ppis_y = binned_ppis.max(axis=1)
    avrgs = gene_property.groupby(binning).median()
    xlim = (n_bins - 0.5, -0.5)  # flip direction of x axis
    if ylabel is None and gene_property.name is not None:
        ylabel = gene_property.name.replace('_', '\n').capitalize()
    panel_ratios = [1 / size_ratio - 1., 1]
    gs = gridspec.GridSpecFromSubplotSpec(2, 2,
                                          subplot_spec=ax.get_subplotspec(),
                                          hspace=pad, wspace=pad,
                                          width_ratios=panel_ratios,
                                          height_ratios=panel_ratios[::-1])
    ax_main = plt.subplot(gs.new_subplotspec((1, 0), rowspan=1, colspan=1))
    new_axes = [ax_main, ]
    if log:
        vmin = 1.
    else:
        vmin = 0.0001
    if log:
        norm = mpl.colors.LogNorm()
    else:
        norm = None
    counts, _binsx, _binsy = np.histogram2d(binned_ppis_x,
                                            binned_ppis_y,
                                            bins=[i - 0.5 for i in range(n_bins + 1)],)
    if correct_diagonal:
        n = genes_ranked.shape[0] / n_bins
        if self_interactions:
            f = n**2 / (n**2 / 2 + n / 2)
        else:
            f = n**2 / (n**2 / 2 - n / 2)
        counts[np.diag_indices(counts.shape[0])] = counts.diagonal() * f
    img = ax_main.imshow(counts.T,
                         origin='lower',
                         cmap=cmap,
                         vmin=vmin,
                         vmax=vmax,
                         norm=norm)
    ax_main.spines['bottom'].set_visible(False)
    ax_main.spines['left'].set_visible(False)
    ax_main.set_xticks([])
    ax_main.set_yticks([])
    ax_main.set_xlim(xlim)
    if label is not None:
        if label_color is None:
            label_color = color
        ax_main.text(n_bins / 2., n_bins / 3., label,
                     color=label_color, ha='right')
    if draw_up:
        ax_up = plt.subplot(gs.new_subplotspec((0, 0), rowspan=1, colspan=1))
        new_axes.append(ax_up)
        ax_up.bar(x=avrgs.index, height=avrgs.values, width=1.0,
                  color=color)
        ax_up.set_xlim(xlim)
        ax_up.set_xticks([])
        ax_up.spines['top'].set_visible(False)
        ax_up.spines['right'].set_visible(False)
        if logy:
            ax_up.set_yscale('log')
        ax_up.set_ylabel(ylabel)
        ax_up.set_ylim(ylim)
        if yticks is not None:
            ax_up.set_yticks(yticks)
        ax_up.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%g'))
    if draw_right:
        ax_rt = plt.subplot(gs.new_subplotspec((1, 1), rowspan=1, colspan=1))
        new_axes.append(ax_rt)
        ax_rt.barh(y=avrgs.index, width=avrgs.values, height=1.0,
                   color=color)
        ax_rt.set_ylim((-0.5, n_bins - 0.5))
        ax_rt.set_yticks([])
        ax_rt.spines['bottom'].set_visible(False)
        ax_rt.spines['right'].set_visible(False)
        if logy:
            ax_rt.set_xscale('log')
        ax_rt.set_xlabel(ylabel)
        ax_rt.xaxis.set_label_position('top')
        ax_rt.set_xlim(ylim)
        if yticks is not None:
            ax_rt.set_xticks(yticks)
        ax_rt.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%g'))
        ax_rt.xaxis.tick_top()
    if draw_scale:
        cax = inset_axes(ax_main,
                         width='{:.0%}'.format(colorbar_width),
                         height='10%',
                         loc='lower left')
        if zticks is None:
            if log:
                zticks = [1, 100, 1000, 10000]
            else:
                zticks = [0, round(counts.max() / 2., 0), counts.max()]
        cbar_label = 'Number of interactions'
        cb = plt.colorbar(img,
                          cax=cax,
                          orientation='horizontal',
                          label=cbar_label,
                          ticks=zticks,
                          format='%d')
        if len(zticks) > 0:
            if zticks[0] == 0:
                cb.set_ticks([vmin] + zticks[1:])
                cb.set_ticklabels([int(zt) for zt in zticks])
    return new_axes, counts


def samogram_double(ppis, gene_property, n_bins,
                    ax=None,
                    id_a=None, id_b=None,
                    correct_diagonal=False,
                    self_interactions=True,
                    draw_up=True, draw_right=False, draw_scale=True,
                    reverse=False,
                    vmax=(None, None),
                    ylim=None, ylabel=None, yticks=None,
                    log=False, logy=False,
                    zticks=(None, None),
                    cmaps=('Purples', 'Blues'),
                    color='grey',
                    size_ratio=0.1,
                    pad=0.04,
                    colorbar_width=0.25,
                    labels=('', '')):
    """2D histogram of PPIs with genes ordered by a given property
    Bar charts of the median of the property per bin can optionally be drawn.
    Gene property values are shuffled before sorting, so that genes with equal
    values are in a random order.
    See the :ref:`tutorial </samogram.ipynb>` for more information.
    Note:
        The diagonal bins have approximately half the number of possible
        interactions as the other bins, since the interactions are undirected.
        Specifically they have :math:`n^2/2 + n/2` possible interactions, with
        the other bins having :math:`n^2`, where :math:`n` is the number of
        genes per bin. This effect is counteracted by the much higher
        likelihood for genes to self-interact than interact with a
        random other gene. The smaller the bin size, the more
        self-interactions will have a visible impact.
    Args:
        ppis (tuple(pandas.DataFrame, pandas.DataFrame)): two tables of interactions, one row per pair
        gene_property (pandas.Series): per-gene quantity to order proteins by
        n_bins (int): number of bins, heatmap will be n_bins x n_bins
        ax (matplotlib.axes.Axes): Axes to draw plot onto
        id_a/id_b (str): name of column containing gene/protein identifiers
        correct_diagonal (bool): increase number of interactions in diagonal
                                 bins to account for there being less possible
                                 combinations. Correction is dependent on value
                                 of `self_interactions` argument.
        self_interactions (bool): whether self interactions are plotted.
        draw_up/draw_right (bool): whether to add bar charts above / to the right
        draw_scale (bool): whether to draw the color scale
        reverse (bool): if true flips direction of ranking of genes
        vmax (int): upper limits on heatmap
        ylim ((float, float)): bar chart axis limits
        ylabel (str): bar chart axis label
        yticks (list(float)): bar chart axis ticks
        zticks (tuple(list(float), list(float))): color scale axis ticks
        log/logy (bool): log scale for heatmap / bar charts
        cmaps (tuple(str, str)): colormap for heatmap
        color (str): color for bar charts
        size_ratio (float): fraction of size of bar charts to heatmap
        pad (float): space between bar charts and heatmap
        colorbar_width (float): size of colorbar as fraction of axes width
        labels (tuple(str, str)): label to put next to heatmap
    Returns:
        (list(matplotlib.axes.Axes)): new axes created
    See Also:
        :func:`samogram`: individual samogram
    Examples:
        Plot samogram of HI-III with default options:
        .. plot::
            :context: close-figs
            >>> import ccsblib.ccsbplotlib as cplt
            >>> from ccsblib import huri
            >>> hi_iii = huri.load_nw_hi_iii(id_type='ensembl_gene_id')
            >>> lit_bm = huri.load_nw_lit_bm_17(id_type='ensembl_gene_id')
            >>> n_pub = huri.load_number_publications_per_gene()
            >>> cplt.samogram_double((hi_iii, lit_bm),
            ...                      n_pub,
            ...                      n_bins=40,
            ...                      id_a='ensembl_gene_id_a',
            ...                      id_b='ensembl_gene_id_b',
            ...                      color='grey',
            ...                      draw_right=True,
            ...                      cmaps=('Purples', 'Blues'),
            ...                      zticks=([0, 100, 200],
            ...                              [0, 100, 200]),
            ...                      vmax=(200, 200),
            ...                      ylabel='Number of publications',
            ...                      labels=('HI-III', 'Lit-BM'))
            >>> import matplotlib.pyplot as plt; fig = plt.gcf(); fig.set_size_inches(6, 6)
    """
    if size_ratio <= 0. or size_ratio >= 1.:
        raise ValueError('size_ratio must be between 0 and 1')
    if ax is None:
        ax = plt.gca()
    if id_a is None:
        id_a = gene_property.index.name + '_a'
    if id_b is None:
        id_b = gene_property.index.name + '_b'
    if isinstance(vmax, (float, int)):
        vmax = (vmax, vmax)
    if isinstance(cmaps, str):
        cmaps = (cmaps, cmaps)
    if isinstance(zticks[0], (float, int)):
        zticks = (zticks, zticks)
    cmap_a = plt.get_cmap(cmaps[0])
    cmap_b = plt.get_cmap(cmaps[1])
    cmap_a.set_under(color=(0, 0, 0, 0))  # fully transparent
    cmap_b.set_under(color=(0, 0, 0, 0))
    genes_ranked = (gene_property.sample(frac=1)  # randomize order first
                                 .sort_values()
                                 .to_frame())
    genes_ranked['rank'] = genes_ranked.reset_index().index.values
    binning = pd.qcut(genes_ranked['rank'], n_bins, labels=False)
    if reverse:
        binning = (n_bins - 1) - binning
    binned_ppis_a = (ppis[0].loc[ppis[0][id_a].isin(gene_property.index) &
                                 ppis[0][id_b].isin(gene_property.index) &
                                 ((ppis[0][id_a] != ppis[0][id_b]) | self_interactions),
                                 [id_a, id_b]]
                            .apply(lambda x: x.map(binning), axis=0))
    binned_ppis_a_x = binned_ppis_a.min(axis=1)
    binned_ppis_a_y = binned_ppis_a.max(axis=1)
    binned_ppis_b = (ppis[1].loc[ppis[1][id_a].isin(gene_property.index) &
                                 ppis[1][id_b].isin(gene_property.index) &
                                 ((ppis[1][id_a] != ppis[1][id_b]) | self_interactions),
                                 [id_a, id_b]]
                            .apply(lambda x: x.map(binning), axis=0))
    binned_ppis_b_x = binned_ppis_b.min(axis=1)
    binned_ppis_b_y = binned_ppis_b.max(axis=1)
    avrgs = gene_property.groupby(binning).median()
    xlim = (n_bins - 0.5, -0.5)  # flip direction of x axis
    if ylabel is None and gene_property.name is not None:
        ylabel = gene_property.name.replace('_', '\n').capitalize()
    panel_ratios = [1 / size_ratio - 1., 1]
    # cut out sqaure sub-area of given Axes
    # get size in inches and convert to figure fraction
    fig = plt.gcf()
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width_scale = min(bbox.width, bbox.height) / bbox.width
    height_scale = min(bbox.width, bbox.height) / bbox.height
    x0, y0, x1, y1 = ax.figbox.get_points().flatten()
    x0 = x0 + (x1 - x0) * (1 - width_scale) / 2
    x1 = x1 - (x1 - x0) * (1 - width_scale) / 2
    y0 = y0 + (y1 - y0) * (1 - height_scale) / 2
    y1 = y1 - (y1 - y0) * (1 - height_scale) / 2
    gs = gridspec.GridSpec(2, 2,
                           left=x0,
                           right=x1,
                           bottom=y0,
                           top=y1,
                           hspace=pad, wspace=pad,
                           width_ratios=panel_ratios,
                           height_ratios=panel_ratios[::-1])
    ax_main = plt.subplot(gs[1, 0])
    new_axes = [ax_main]
    if log:
        vmin = 1.
    else:
        vmin = 0.0001
    if log:
        norm = mpl.colors.LogNorm()
    else:
        norm = None
    counts_a, _binsx, _binsy = np.histogram2d(binned_ppis_a_x,
                                              binned_ppis_a_y,
                                              bins=[i - 0.5 for i in range(n_bins + 1)],)
    # Note the reversal of x and y below
    counts_b, _binsx, _binsy = np.histogram2d(binned_ppis_b_y,
                                              binned_ppis_b_x,
                                              bins=[i - 0.5 for i in range(n_bins + 1)],)
    if correct_diagonal:
        n = genes_ranked.shape[0] / n_bins
        if self_interactions:
            f = n**2 / (n**2 / 2 + n / 2)
        else:
            f = n**2 / (n**2 / 2 - n / 2)
        counts_a[np.diag_indices(counts_a.shape[0])] = counts_a.diagonal() * f
        counts_b[np.diag_indices(counts_b.shape[0])] = counts_b.diagonal() * f
    img_a = ax_main.imshow(counts_a.T,
                           origin='lower',
                           cmap=cmap_a,
                           vmin=vmin,
                           vmax=vmax[0],
                           norm=norm)
    img_b = ax_main.imshow(counts_b.T,
                           origin='lower',
                           cmap=cmap_b,
                           vmin=vmin,
                           vmax=vmax[1],
                           norm=norm)
    xhigh, xlow = ax_main.get_xlim()
    ylow, yhigh = ax_main.get_ylim()
    crop_triangle_a = mpl.patches.Polygon([[xlow, yhigh],
                                           [xhigh, ylow],
                                           [xhigh, yhigh]],
                                          transform=ax_main.transData)
    img_a.set_clip_path(crop_triangle_a)
    crop_triangle_b = mpl.patches.Polygon([[xlow, yhigh],
                                           [xhigh, ylow],
                                           [xlow, ylow]],
                                          transform=ax_main.transData)
    img_b.set_clip_path(crop_triangle_b)
    ax_main.set_xticks([])
    ax_main.set_yticks([])
    ax_main.set_xlim(xlim)
    # diagonal dividing line
    ax_main.plot(ax_main.get_xlim()[::-1],
                 ax_main.get_ylim(),
                 color=ax_main.spines['top'].get_edgecolor(),
                 linewidth=ax_main.spines['top'].get_linewidth(),
                 linestyle=ax_main.spines['top'].get_linestyle())
    if draw_up:
        ax_up = plt.subplot(gs[0, 0], sharex=ax_main)
        new_axes.append(ax_up)
        ax_up.bar(x=avrgs.index, height=avrgs.values, width=1.0,
                  color=color)
        ax_up.set_xlim(xlim)
        ax_up.set_xticks([])
        ax_up.spines['top'].set_visible(False)
        ax_up.spines['right'].set_visible(False)
        if logy:
            ax_up.set_yscale('log')
        ax_up.set_ylabel(ylabel)
        ax_up.set_ylim(ylim)
        if yticks is not None:
            ax_up.set_yticks(yticks)
        ax_up.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%g'))
    if draw_right:
        ax_rt = plt.subplot(gs[1, 1], sharey=ax_main)
        new_axes.append(ax_rt)
        ax_rt.barh(y=avrgs.index, width=avrgs.values, height=1.0,
                   color=color)
        ax_rt.set_ylim((-0.5, n_bins - 0.5))
        ax_rt.set_yticks([])
        ax_rt.spines['bottom'].set_visible(False)
        ax_rt.spines['right'].set_visible(False)
        if logy:
            ax_rt.set_xscale('log')
        ax_rt.set_xlabel(ylabel)
        ax_rt.xaxis.set_label_position('top')
        ax_rt.set_xlim(ylim)
        if yticks is not None:
            ax_rt.set_xticks(yticks)
        ax_rt.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%g'))
        ax_rt.xaxis.tick_top()
    if draw_scale:
        zticks_a, zticks_b = zticks
        same_scale = vmax[0] == vmax[1]
        # I don't understand why I have to put 0.25 instead of 0.1....
        cax_a = inset_axes(ax_main,
                           width='{:.0%}'.format(colorbar_width),
                           height='10%',
                           loc='lower center',
                           bbox_to_anchor=[0., -0.35, 1., 1.],
                           bbox_transform=ax_main.transAxes)
        cax_b = inset_axes(ax_main,
                           width='{:.0%}'.format(colorbar_width),
                           height='10%',
                           loc='lower center',
                           bbox_to_anchor=[0., -0.25, 1., 1.],
                           bbox_transform=ax_main.transAxes)
        if zticks_a is None:
            if log:
                zticks_a = [1, 100, 1000, 10000]
            else:
                zticks_a = [0, round(counts_a.max() / 2., 0), counts_a.max()]
        if zticks_b is None:
            if log:
                zticks_b = [1, 100, 1000, 10000]
            else:
                zticks_b = [0, round(counts_b.max() / 2., 0), counts_b.max()]
        if same_scale:
            zticks_b = []
        cb_a = plt.colorbar(img_a,
                            cax=cax_a,
                            orientation='horizontal',
                            label='Number of interactions',
                            ticks=zticks_a,
                            format='%d')
        cb_b = plt.colorbar(img_b,
                            cax=cax_b,
                            orientation='horizontal',
                            label='',
                            ticks=zticks_b)
        cax_b.xaxis.set_ticks_position('top')
        cax_a.text(cb_a.get_clim()[1] * 1.05, 1, labels[0])
        cax_b.text(cb_b.get_clim()[1] * 1.05, 1, labels[1])
        if len(zticks_a) > 0:
            if zticks_a[0] == 0:
                cb_a.set_ticks([vmin] + zticks_a[1:])
            tick_labels = [str(int(zt)) for zt in zticks_a]
            if vmax[0] is not None:
                if vmax[0] < counts_a.max() or ((vmax[0] < counts_b.max()) and same_scale):
                    tick_labels = tick_labels[:-1] + ['≥' + tick_labels[-1]]
            cb_a.set_ticklabels(tick_labels)
        if len(zticks_b) > 0 and not same_scale:
            if zticks_b[0] == 0:
                cb_b.set_ticks([vmin] + zticks_b[1:])
            tick_labels = [str(int(zt)) for zt in zticks_b]
            if vmax[1] is not None:
                if vmax[1] < counts_b.max() and not same_scale:
                    tick_labels = tick_labels[:-1] + ['≥' + tick_labels[-1]]
            cb_b.set_ticklabels(tick_labels)
    return new_axes


def degree_distribution_plot(nw,
                             id_a, id_b,
                             ax=None,
                             ymin=None, xmax=None,
                             **kwargs):
    """Degree distribution on log-log scale with log binning.
    Based on the plots and advice in http://barabasi.com/networksciencebook
    See: Chapter 4, Advanced Topic 3.B
    See the :ref:`tutorial </degree_plot.ipynb>` for more information.
    Args:
        nw (DataFrame): PPI dataset, one row per pair
        id_a/id_b: column names of protein identifiers
        ax (matplotlib.axes.Axes): Axes to draw plot onto
        ymin/xmax (float): axis limits
        **kwargs: passed to ax.plot
    Examples:
        Plot degree distribution of HI-III:
        .. plot::
            :context: close-figs
            >>> import ccsblib.ccsbplotlib as cplt
            >>> from ccsblib import huri
            >>> hi_iii = huri.load_nw_hi_iii(id_type='ensembl_gene_id')
            >>> cplt.degree_distribution_plot(hi_iii,
            ...                               id_a='ensembl_gene_id_a',
            ...                               id_b='ensembl_gene_id_b')
    """
    if ax is None:
        ax = plt.gca()
    d = degree_per_protein(nw, id_a, id_b)
    bins = []
    i = 0
    while 2**i - 0.5 < d.max():
        bins.append(2**i - 0.5)
        i += 1
    pk = ((d.groupby(pd.cut(d, bins)).size().values /
           np.array([bins[j + 1] - bins[j] for j in range(len(bins) - 1)])) /
          d.shape[0])
    dMeans = d.groupby(pd.cut(d, bins)).mean().values
    ax.plot(dMeans, pk, 'o', **kwargs)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim((ymin, 1.0))
    ax.set_xlim((0.8, xmax))
    ax.set_ylabel('Probability density')
    ax.set_xlabel('Degree')


def checkerboard(data,
                 protein_a_column=None, protein_b_column=None,
                 detection_columns=None,
                 sort=True,
                 assay_labels=None,
                 positive_color='yellow', negative_color='white',
                 ax=None):
    """Plot yes/no detection for benchmark PPI set with different assays
    See Braun et al, Nature Methods, 2010 for examples.
    Args:
        data (pandas.DataFrame): PPI test results. No missing values.
        protein_a/b_column (str): columns with protein names
        detection_columns (list(str)): name of columns containing boolean results
        sort (bool): whether to sort pairs by number of assays detected and assay order
        assay_labels (list(str)): names of assays to print
        positive_color (str/RGB/RGBA or list(colors)): single color or list of colors for each different assay
        negative_color (str/RGB/RGBA): color to indicate undetected pairs
        ax (matplotlib.axes.Axes): Axes to draw plot onto
    Examples:
        Make a checkerboard of some dummy data:
        .. plot::
            :context: close-figs
            >>> import pandas as pd
            >>> import ccsblib.ccsbplotlib as cplt
            >>> prs_results = pd.DataFrame(columns=['gene_a', 'gene_b', 'Y2H', 'MAPPIT', 'GPCA'],
            ...                            data=[['ABC1', 'ABC2', False, False, False],
            ...                                  ['EFG1', 'EFG2', False, False, False],
            ...                                  ['HIJ1', 'HIJ2', True, False, False],
            ...                                  ['KLM1', 'KLM2', False, False, False],
            ...                                  ['NOP1', 'NOP2', True, True, True],
            ...                                  ['QRS1', 'QRS2', True, False, True],
            ...                                  ['TUV1', 'TUV2', False, False, True],
            ...                                  ['XYZ1', 'XYZ2', False, False, False]])
            >>> cplt.checkerboard(data=prs_results,
            ...                   protein_a_column='gene_a',
            ...                   protein_b_column='gene_b',
            ...                   detection_columns=['Y2H',
            ...                                      'MAPPIT',
            ...                                      'GPCA'])
    """
    df = data.copy()
    if ax is None:
        ax = plt.gca()
    if assay_labels is None:
        assay_labels = detection_columns
    if protein_a_column is None:
        protein_a_column = df.columns[0]
    if protein_b_column is None:
        protein_a_column = df.columns[1]
    if detection_columns is None:
        detection_columns = list(df.columns[df.dtypes == bool])
    elif isinstance(detection_columns, str):
        detection_columns = [detection_columns]
    df['total_positives'] = df[detection_columns].sum(axis=1)
    if sort:
        df = df.sort_values(by=['total_positives'] + detection_columns, ascending=False)
    if isinstance(positive_color, list) and len(positive_color) == len(detection_columns):
        results = df[detection_columns].values.T
        for i, color in enumerate(positive_color):
            m = np.zeros(shape=results.shape, dtype=bool)
            m[i, :] = True
            ax.imshow(results * m,
                      cmap=mpl.colors.ListedColormap([negative_color if i == 0 else (0, 0, 0, 0),
                                                             color]))
    else:
        ax.imshow(df[detection_columns].values.T,
                  cmap=mpl.colors.ListedColormap([negative_color, positive_color]))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks(range(len(detection_columns)))
    ax.set_yticklabels(assay_labels)
    ax.yaxis.set_tick_params(length=0)
    ax.set_xticks([])
    len_longest_name = df[protein_a_column].str.len().max()
    for i, (name_a, name_b) in enumerate(zip(df[protein_a_column].values, df[protein_b_column].values)):
        ax.text(i, -0.6, name_a + ' ' * (len_longest_name - len(name_a) + 2) + name_b,
                rotation=90,
                fontfamily='monospace',
                va='bottom',
                ha='center')