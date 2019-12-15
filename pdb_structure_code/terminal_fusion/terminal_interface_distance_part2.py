import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
pymol.finish_launching()
from pymol import cmd, stored
from glob import glob
from time import gmtime, strftime

import os
import lib_aux
from lib_aux import *

# paths and output folders
PATH_PDB_PPI_MODELS            = PATH_TO_INTERACTOME3D_DOWNLOADED_INTERACTION_MODELS
PATH_PDB_PPI_EXP               = PATH_TO_INTERACTOME3D_DOWNLOADED_INTERACTION_EXPERIMENTAL_STRUCTURES
PATH_PDB_PPI_DDM               = PATH_TO_INTERACTOME3D_DOWNLOADED_INTERACTION_DOMAIN_BASED_MODELS
OUTPUT_FOLDER_PNG              = "png_files.nosync"
OUTPUT_FOLDER_SESSION          = "pymol_sessions.nosync"

if not os.path.exists(OUTPUT_FOLDER_PNG):
    os.makedirs(OUTPUT_FOLDER_PNG)

if not os.path.exists(OUTPUT_FOLDER_SESSION):
    os.makedirs(OUTPUT_FOLDER_SESSION)

# get all pdb files in folder
cmd.set("dot_solvent", "on")

# open output files
output_distances    = open("output_distances_monomer.csv", "w")
output_ignored      = open("output_ignored_monomer.txt", "w")

# headers of output files
output_distances.write("complex_PDB\t" \
                       "mo_d_N1\tmo_d_C1\tmo_d_N2\tmo_d_C2\n")

input_file = "input_distances_monomer.csv"
input_data = open(input_file, "r")
input_data.readline() # skip first line
lines = input_data.readlines()

# interface residues mapped on monomers
input_iface_monomer_file = "iface_residues_on_monomers.txt"
input_iface_monomer_data = open(input_iface_monomer_file, "r").readlines()
iface_res_monomers = {}
for line in input_iface_monomer_data:
  (key, val) = line.strip().split("\t")
  iface_res_monomers[key] = val.split(";")

# folder to put temp PDB files
tmp_folder = "tmp_pdb.nosync"
if not os.path.exists(tmp_folder):
    os.makedirs(tmp_folder)

done = 0
ignored = 0
i_line = 0

for line in lines:
  i_line += 1
  this_line       = line.strip().split("\t")
  this_co_PDB     = this_line[0]  # complex filename
  this_mo_N1_PDB   = this_line[1] # monomer PDB for Nterminal 1
  this_mo_N1_type  = this_line[2] # monomer PDB for Nterminal 1 (type)
  this_mo_N1_iface = this_line[3] # monomer PDB interface residue
  this_mo_N1_res   = this_line[4] # monomer PDB Nterminal residue (terminal residue)
  this_mo_C1_PDB   = this_line[5] # same for Cterminal 1
  this_mo_C1_type  = this_line[6]
  this_mo_C1_iface = this_line[7]
  this_mo_C1_res   = this_line[8]
  this_mo_N2_PDB   = this_line[9] # same for Nterminal 2
  this_mo_N2_type  = this_line[10]
  this_mo_N2_iface = this_line[11] 
  this_mo_N2_res   = this_line[12]
  this_mo_C2_PDB   = this_line[13] # same for Cterminal 2
  this_mo_C2_type  = this_line[14]
  this_mo_C2_iface = this_line[15] 
  this_mo_C2_res   = this_line[16]

  try:
      mo_N1_iface_residues = iface_res_monomers[ "%s_N1" % this_co_PDB ]
      mo_C1_iface_residues = iface_res_monomers[ "%s_C1" % this_co_PDB ]
      mo_N2_iface_residues = iface_res_monomers[ "%s_N2" % this_co_PDB ]
      mo_C2_iface_residues = iface_res_monomers[ "%s_C2" % this_co_PDB ]
      ### MONOMER distances
      mo_N1 = f_get_monomer_distances(this_mo_N1_PDB, this_mo_N1_type, this_mo_N1_res, this_mo_N1_iface, mo_N1_iface_residues)
      mo_C1 = f_get_monomer_distances(this_mo_C1_PDB, this_mo_C1_type, this_mo_C1_res, this_mo_C1_iface, mo_C1_iface_residues)
      mo_N2 = f_get_monomer_distances(this_mo_N2_PDB, this_mo_N2_type, this_mo_N2_res, this_mo_N2_iface, mo_N2_iface_residues)
      mo_C2 = f_get_monomer_distances(this_mo_C2_PDB, this_mo_C2_type, this_mo_C2_res, this_mo_C2_iface, mo_C2_iface_residues)

      output_distances.write("%s" % (this_co_PDB))
      output_distances.write("\t%.1f"*4 % (float(mo_N1['e_dist']), float(mo_C1['e_dist']), float(mo_N2['e_dist']), float(mo_C2['e_dist']) ) )
      output_distances.write("\n")
      done += 1

  except:
    print "%s: exception\n" % this_co_PDB
    ignored += 1
    output_ignored.write("%s\n" % this_co_PDB)

  cmd.delete("*")

  if (i_line % 10 == 0):
   print "Done %s/%s (%s)\n" % (i_line, len(lines), strftime("%Y-%m-%d %H:%M:%S", gmtime()))

output_distances.close()
output_ignored.close()

print done


