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
output_distances = open("output_distances_complex.csv", "w")
output_ignored   = open("output_ignored_PDBs_complex.txt", "w") # PDBs not analyzed (not enough atoms different than CA)
output_int_res   = open("output_interface_residues.txt", "w")   # output interface residue data

# headers of output files
output_distances.write("complex_PDB\tPDB_type\tiface_res1\tiface_dist1\tiface_res2\tiface_dist2\t" \
                       "co_d_N1\tco_d_C1\tco_d_N2\tco_d_C2\n") #t"     \

input_file = "input_distances_complex.csv"
input_data = open(input_file, "r")
input_data.readline() # skip first line
lines = input_data.readlines()

# folder to put temp PDB files
tmp_folder = "tmp_pdb.nosync"
if not os.path.exists(tmp_folder):
    os.makedirs(tmp_folder)

done = 0
ignored = 0
i_line = 0

all_iface_data = [] # to store the interface residues for each chain of the complex
for line in lines:
  i_line += 1
  this_line       = line.strip().split("\t")
  this_co_PDB     = this_line[0]  # complex filename
  this_co_type    = this_line[1]  # complex type of PDB: model or structure (important to get the data folder w/ the PDBs)
  this_co_N_res1  = this_line[2]  # complex N_term residue of protein1
  this_co_C_res1  = this_line[3]  # complex C_term residue of protein1
  this_co_N_res2  = this_line[4]  # complex N_term residue of protein2
  this_co_C_res2  = this_line[5]  # complex C_term residue of protein2

  if (this_co_type == "Model"):
    this_co_PDB_path = PATH_PDB_PPI_MODELS
  elif (this_co_type == "Structure"):
    this_co_PDB_path = PATH_PDB_PPI_EXP
  else:
    this_co_PDB_path = PATH_PDB_PPI_DDM
    

  this_PDB_file = "%s/%s" % (this_co_PDB_path, this_co_PDB)

  try:
    # load pdb and calculate interface on protein A and B (only CA atoms)
    cmd.load(this_PDB_file, "this_complex")

    # we require the PDB to have at least 10 atoms different than CA
    if (len(cmd.get_model("not name ca").atom) > 10):
      iface_data = f_get_interface_center("this_complex")
      all_iface_data.append([this_co_PDB, "A", ";".join(iface_data['int_res1']), "B", ";".join(iface_data['int_res2'])] )

      ### distance from complex terminals to interface centers -> we use "first"
      co_d_N_iface1 = cmd.get_distance("first this_complex & chain A & name ca & resi %s" % this_co_N_res1, "first this_complex & chain A & name ca & resi %s" % iface_data['CA1'])
      co_d_C_iface1 = cmd.get_distance("first this_complex & chain A & name ca & resi %s" % this_co_C_res1, "first this_complex & chain A & name ca & resi %s" % iface_data['CA1'])
      co_d_N_iface2 = cmd.get_distance("first this_complex & chain B & name ca & resi %s" % this_co_N_res2, "first this_complex & chain B & name ca & resi %s" % iface_data['CA2'])
      co_d_C_iface2 = cmd.get_distance("first this_complex & chain B & name ca & resi %s" % this_co_C_res2, "first this_complex & chain B & name ca & resi %s" % iface_data['CA2'])

      output_distances.write("%s\t%s\t%s\t%.1f\t%s\t%.1f" % (this_co_PDB, this_co_type, iface_data['CA1'], iface_data['d_CA1'], iface_data['CA2'], iface_data['d_CA2']) )
      output_distances.write("\t%.1f"*4 % (float(co_d_N_iface1), float(co_d_C_iface1), float(co_d_N_iface2), float(co_d_C_iface2) ))
      output_distances.write("\n")
      done += 1

    else:
      print "\n%s: ignored\n" % this_co_PDB
      output_ignored.write("%s\n" % this_co_PDB)

  except:
    print "%s: exception\n" % this_co_PDB
    ignored += 1
    output_ignored.write("%s\n" % this_co_PDB)

  cmd.delete("*")

  if (i_line % 10 == 0):
   print "Done %s/%s (%s)\n" % (i_line, len(lines), strftime("%Y-%m-%d %H:%M:%S", gmtime()))

output_distances.close()
output_ignored.close()

# output interface residue data
for item in all_iface_data:
  output_int_res.write("%s\n" % item)

print done


