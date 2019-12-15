import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
pymol.finish_launching()
from pymol import cmd, stored
from glob import glob


interaction_file = "interactions.dat"
interaction_data = open(interaction_file, "r").read().split("\n")
pdb_path = {"Structure": PDB_EXPERIMENTAL_PATH,
     	    "Model": PDB_MODEL_PATH,
            "Dom_dom_model": PDB_DOMAIN_BASED_PATH }

cmd.set("dot_solvent", "on")
output = open("interface_area.csv", "w")
output.write("protein1\tprotein2\tPDB\tPDB_type\tiface_area\tcomplex_area\tP1_area\tP2_area\n")
for this_line in interaction_data[1:]:
 if (this_line != ""):
  this_info = this_line.split("\t")
  this_p1   = this_info[0]
  this_p2   = this_info[1]
  this_type = this_info[4]
  this_pdb  = this_info[21]
  this_file = "%s%s" % (pdb_path[this_type], this_pdb)
  print this_file
  cmd.load(this_file, "this_pdb")
  cmd.copy("ca", "this_pdb")
  cmd.remove("ca and chain B")
  cmd.copy("cb", "this_pdb")
  cmd.remove("cb and chain A")
  complex_area = cmd.get_area("this_pdb")
  chainA_area = cmd.get_area("ca")
  chainB_area = cmd.get_area("cb")
  interface_area = (chainA_area + chainB_area - complex_area)/2
  output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (this_p1, this_p2, this_pdb, this_type, interface_area, complex_area, chainA_area, chainB_area))
  cmd.delete("this_pdb")
  cmd.delete("ca")
  cmd.delete("cb")
output.close()

