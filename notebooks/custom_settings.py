"""
Holds color and label definitions and matplotlib settings.
"""
HURI_COLOR = (155 / 255, 97 / 255, 153 / 255)
LIT_COLOR = (60 / 255, 134 / 255, 184 / 255)
# from colorbrewer; middle three of 5-class purples
ASSAY_V1_COLOR = '#cbc9e2'
ASSAY_V2_COLOR = '#9e9ac8'
ASSAY_V3_COLOR = '#756bb1'
HI_I_05_COLOR = 'pink'
HI_II_14_COLOR = 'orchid'
RRS_COLOR = (224 / 255, 16 / 255, 28 / 255)
# from colorbrewer; middle three of 5-class oranges
BIOPLEX_COLOR = '#e6550d'
QUBIC_COLOR = '#fd8d3c'
COFRAC_COLOR = '#fdbe85'
FORMATS = ['.pdf', '.png']
HI_I_to_III_TICK_LABELS = ['HI-I-05', 'HI-II-14', 'HI-III-19']
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
hi_cmap = LinearSegmentedColormap.from_list('HI-III-19', ['white', HURI_COLOR, 'black'], 1000)
lit_cmap = LinearSegmentedColormap.from_list('Lit-BM', ['white', LIT_COLOR, 'black'], 1000)
plt.register_cmap(name=hi_cmap.name, cmap=hi_cmap)
plt.register_cmap(name=lit_cmap.name, cmap=lit_cmap)

"""
Path and database-related variables
"""
#DB_NAME = 'hi2018_paper'
DB_NAME = 'bioinfo_katja'

"""
Used cutoffs
"""
GTEX_EXPR_CUTOFF = 5
TIP_CUTOFF = 2

"""
network names
"""
NETWORK_NAMES = ['HI-III','HI-II-14','HI-I-05','Lit-BM-17','BioPlex','QUBIC','CoFrac']
