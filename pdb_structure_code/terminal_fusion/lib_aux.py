import pymol
from pymol import cmd, stored
import interfaceResidues
import center_of_mass
import distancetoatom
import os

PATH_PDB_MONO_MOD              = PATH_TO_INTERACTOME3D_HURI_RESULT_PROTEIN_MODEL_PDBS
PATH_PDB_MONO_EXP              = PATH_TO_INTERACTOME3D_HURI_RESULT_PROTEIN_EXPERIMENTAL_PDBS

if not os.path.exists(OUTPUT_FOLDER_SURFACE_DISTANCE):
    os.makedirs(OUTPUT_FOLDER_SURFACE_DISTANCE)

def output_euclidian_distances_header(output_euclidian):
  output_euclidian.write("PDB_file\t" \
                         "d_AN_to_A_iface\td_AN1_to_A_iface\t" \
                         "d_AC_to_A_iface\td_AC1_to_A_iface\t" \
                         "d_BN_to_B_iface\td_BN1_to_B_iface\t" \
                         "d_BC_to_B_iface\td_BC1_to_B_iface\t" \
                         "d_AN_BN\td_AN_BC\td_AC_BN\td_AC_BC\td_AN_AC\td_BN_BC\n")


def output_euclidian_distances(output_euclidian, this_PDB_file, \
                               euc_dist_Aterms_iface_0, euc_dist_Aterms_iface_1, \
                               euc_dist_Bterms_iface_0, euc_dist_Bterms_iface_1, \
                               euc_dist_bw_terminals):
  output_euclidian.write("%s" % this_PDB_file)
  output_euclidian.write("\t%.1f"*14 %
                         (euc_dist_Aterms_iface_0['d_N_term_iface'], euc_dist_Aterms_iface_1['d_N_term_iface'], \
                          euc_dist_Aterms_iface_0['d_C_term_iface'], euc_dist_Aterms_iface_1['d_C_term_iface'], \
                          euc_dist_Bterms_iface_0['d_N_term_iface'], euc_dist_Bterms_iface_1['d_N_term_iface'], \
                          euc_dist_Bterms_iface_0['d_C_term_iface'], euc_dist_Bterms_iface_1['d_C_term_iface'], \
                          euc_dist_bw_terminals['d_A_N_term_to_B_N_term'], euc_dist_bw_terminals['d_A_N_term_to_B_C_term'], \
                          euc_dist_bw_terminals['d_A_C_term_to_B_N_term'], euc_dist_bw_terminals['d_A_C_term_to_B_C_term'], \
                          euc_dist_bw_terminals['d_A_N_term_to_A_C_term'], euc_dist_bw_terminals['d_B_N_term_to_B_C_term']) )
  output_euclidian.write("\n")


def output_protein_info_header(output_protein_info):
  output_protein_info.write("PDB_file\t" \
                            "A_n_res\tA_n_res_iface\tA_center_iface\tA_center_iface_offset\tA_center_iface_far\tA_center_iface_far_dist\t" \
                            "A_N_term\tA_N_term1\tA_C_term\tA_C_term1\t" \
                            "B_n_res\tB_n_res_iface\tB_center_iface\tB_center_iface_offset\tB_center_iface_far\tB_center_iface_far_dist\t" \
                            "B_N_term\tB_N_term1\tB_C_term\tB_C_term1\n")


# protein info: n_residues, n_residues_interface, closest_CA_to_interface, farthest_CA_to_interface
def output_protein_info_data(output_protein_info, this_PDB_file, \
                             distances_from_A_interface_center, euc_dist_Aterms_iface_0, euc_dist_Aterms_iface_1, \
                             distances_from_B_interface_center, euc_dist_Bterms_iface_0, euc_dist_Bterms_iface_1):
  # data for chain A
  A_n_res       = cmd.count_atoms("chain A & name CA")
  A_n_res_iface = cmd.count_atoms("A_interface and name CA")
  A_center_iface_CA            = distances_from_A_interface_center[2][0].split("/")[4].split("`")[1]  # get the closest CA
  A_center_farthest_iface_CA   = distances_from_A_interface_center[-1][0].split("/")[4].split("`")[1] # get farthest interface CA residue from the interface center (and distance)
  d_A_center_iface_CA          = "%.1f" % distances_from_A_interface_center[2][4]                     # get the distance
  d_A_center_farthest_iface_CA = "%.1f" % distances_from_A_interface_center[-1][4]
  # data for chain B
  B_n_res       = cmd.count_atoms("chain B & name ca")
  B_n_res_iface = cmd.count_atoms("B_interface and name CA")
  B_center_iface_CA            = distances_from_B_interface_center[2][0].split("/")[4].split("`")[1]  # get the closest CA
  B_center_farthest_iface_CA   = distances_from_B_interface_center[-1][0].split("/")[4].split("`")[1] # get farthest interface CA residue from the interface center (and distance)
  d_B_center_iface_CA          = "%.1f" % distances_from_B_interface_center[2][4]                     # get the distance
  d_B_center_farthest_iface_CA = "%.1f" % distances_from_B_interface_center[-1][4]
  output_protein_info.write("%s" % this_PDB_file)
  output_protein_info.write("\t%s"*20 % (A_n_res, A_n_res_iface, A_center_iface_CA, d_A_center_iface_CA, A_center_farthest_iface_CA, d_A_center_farthest_iface_CA, \
                                         euc_dist_Aterms_iface_0["N_term_res"], euc_dist_Aterms_iface_1["N_term_res"], \
                                         euc_dist_Aterms_iface_0["C_term_res"], euc_dist_Aterms_iface_1["C_term_res"], \
                                         B_n_res, B_n_res_iface, B_center_iface_CA, d_B_center_iface_CA, B_center_farthest_iface_CA, d_B_center_farthest_iface_CA, \
                                         euc_dist_Bterms_iface_0["N_term_res"], euc_dist_Bterms_iface_1["N_term_res"], \
                                         euc_dist_Bterms_iface_0["C_term_res"], euc_dist_Bterms_iface_1["C_term_res"]) )
  output_protein_info.write("\n")


# iface_res: central interface residue
# iface_residues: interface residues
def f_get_monomer_distances(monomer_PDB, monomer_PDB_type, term_residue, iface_res, iface_residues):
  e_dist = -1
  s_dist = -1
  min_iface_dist = -1
  try:
   if ((monomer_PDB != "NA") & (iface_res != "NA") & (term_residue != "NA")):
    if (monomer_PDB_type == "Model"):
      this_PDB_path = PATH_PDB_MONO_MOD
    else:
      this_PDB_path = PATH_PDB_MONO_EXP
    PDB_file = "%s/%s" % (this_PDB_path, monomer_PDB)
    cmd.load(PDB_file, "this_mono")
    # distance from terminal to interface center
    e_dist = cmd.get_distance("first this_mono & name ca & resi %s" % term_residue, "first this_mono & name ca & resi %s" % iface_res)
    ## minimum distance to interface residues
    if (iface_residues != ["NA"]):
      if (term_residue in iface_residues):
        min_iface_dist = 0
      else:
        resi_txt = " or ".join(["resi " + this_res for this_res in iface_residues])
        min_iface_dist = distancetoatom.distancetoatom(origin="first this_mono & name ca & resi %s" % term_residue, selection="this_mono & name ca & (%s)" % resi_txt, property_name="", cutoff=1000)[2][4]
  except:
    print "ERROR in monomer distances %s (%s) %s %s\n" % (monomer_PDB, monomer_PDB_type, term_residue, iface_res)
  cmd.delete("this_mono")
  to_return = { 'e_dist': e_dist }
  return to_return


# Identifies the interfaces of chain A and chain B of complex_txt (named "A_interface" and "B_interface")
# calculates the center of coordinates of interfaces 1 and 2, which are created as pseudoatoms named "A_center_pseudoatom" and "B_center_pseudoatom"
# it also retrieves (and returns) the closest CA to each of those interface centers
def f_get_interface_center(complex_txt):
  chains = cmd.get_chains(complex_txt)
  interface_data = interfaceResidues.interfaceResidues(complex_txt, cA=chains[0], cB=chains[1], cutoff=1, selName="Interface_res")
  cmd.select("A_interface", "Interface_res & chain A")
  cmd.select("B_interface", "Interface_res & chain B")

  # get residues
  stored.intA = []
  cmd.iterate("A_interface & name CA", "stored.intA.append(resi)")
  stored.intB = []
  cmd.iterate("B_interface & name CA", "stored.intB.append(resi)")

  # get center of coordinates of the interface (A and B)
  center_of_mass.com("A_interface", object="A_center_pseudoatom")
  center_of_mass.com("B_interface", object="B_center_pseudoatom")
  A_center = center_of_mass.get_com("A_interface")
  B_center = center_of_mass.get_com("B_interface")
 
  dist_iface_center1 = distancetoatom.distancetoatom(origin="A_center_pseudoatom", selection="A_interface and name CA", property_name="", cutoff=200)
  center_iface_CA1   = dist_iface_center1[2][0].split("/")[4].split("`")[1]  # get the closest CA
  d_center_iface_CA1 = dist_iface_center1[2][4]                              # get the distance

  dist_iface_center2 = distancetoatom.distancetoatom(origin="B_center_pseudoatom", selection="B_interface and name CA", property_name="", cutoff=200)
  center_iface_CA2   = dist_iface_center2[2][0].split("/")[4].split("`")[1]  # get the closest CA
  d_center_iface_CA2 = dist_iface_center2[2][4]                              # get the distance

  to_return = {"CA1"     : center_iface_CA1,
               "d_CA1"   : d_center_iface_CA1,
               "CA2"     : center_iface_CA2,
               "d_CA2"   : d_center_iface_CA2,
               "int_res1": stored.intA,
               "int_res2": stored.intB}
  return to_return


def output_PNG(this_PDB_file, OUTPUT_FOLDER_PNG, OUTPUT_FOLDER_SESSION):
  # show complex
  cmd.hide("all")
  cmd.bg_color("white")
  cmd.show("cartoon", "this_complex")
  # show protein A
  cmd.set_color("A_interface_color", [167, 244, 66]) # lightgreen
  cmd.set_color("A_C_term_color", [252, 245, 146])   # pale yellow
  cmd.color("yellow", "chain A")
  cmd.color("A_interface_color", "A_interface")
  cmd.alter("A_center_pseudoatom", "vdw=2")
  cmd.show("sphere", "A_center_pseudoatom")
  cmd.color("green", "A_center_pseudoatom")
  cmd.show("sphere", "A_N_term")
  cmd.color("yellow", "A_N_term")
  cmd.show("sphere", "A_C_term")
  cmd.color("A_C_term_color", "A_C_term")
  # show protein B
  cmd.set_color("B_interface_color", [80, 181, 244]) # medium blue
  cmd.set_color("B_C_term_color", [179, 229, 229])   # pale cyan
  cmd.color("cyan", "chain B")
  cmd.color("B_interface_color", "B_interface")
  cmd.alter("B_center_pseudoatom", "vdw=2")
  cmd.show("sphere", "B_center_pseudoatom")
  cmd.color("blue", "B_center_pseudoatom")
  cmd.show("sphere", "B_N_term")
  cmd.color("cyan", "B_N_term")
  cmd.show("sphere", "B_C_term")
  cmd.color("B_C_term_color", "B_C_term")
  # nice output
  cmd.rebuild()
  cmd.zoom("all")
  cmd.png("%s/%s.png" % (OUTPUT_FOLDER_PNG, this_PDB_file), quiet=1)
  # save session
  cmd.save("%s/%s.pse" % (OUTPUT_FOLDER_SESSION, this_PDB_file))


## calculates euclidian distance from N-term and C_term of specified chain to "residue"
# position = number of residues away from the terminal
# chain: chain whose terminals have to be calculated
def get_terminal_residues_and_euclidian_distances_from_residue(chain, residue, position=0):
  N_res = len(cmd.get_model("chain %s & name ca" % chain).atom)
  if (position < N_res):
    N_term_res         = cmd.get_model("chain %s & name ca and resi 1-" % chain).atom[position].resi            # smallest residue is 1 (negative residues are ignored)
    C_term_res         = cmd.get_model("chain %s & name ca" % chain).atom[-(position+1)].resi
    d_N_term_iface     = cmd.get_distance("first chain %s & name ca & resi %s" % (chain, N_term_res), residue)  # we have to use "first" because there was some strange error message
    d_C_term_iface     = cmd.get_distance("first chain %s & name ca & resi %s" % (chain, C_term_res), residue)
  else:
    N_term_res = "NA"
    C_term_res = "NA"
    d_N_term_iface = float('nan')
    d_C_term_iface = float('nan')
  to_return = {'N_term_res' : N_term_res,
               'C_term_res' : C_term_res,
               'd_N_term_iface' : d_N_term_iface,
               'd_C_term_iface' : d_C_term_iface}
  return to_return


## selects A and B N- and C- terminals and calculates all distances between them
def get_euclidian_distances_bw_terminals():
  # get N_terms and C_terms (only CA atoms)
  cmd.select("A_N_term", "first chain A & name ca")
  cmd.select("A_C_term", "last chain A & name ca")
  cmd.select("B_N_term", "first chain B & name ca")
  cmd.select("B_C_term", "last chain B & name ca")

  # euclidian distance between terminals
  d_A_N_term_to_B_N_term = cmd.get_distance("A_N_term", "B_N_term")
  d_A_N_term_to_B_C_term = cmd.get_distance("A_N_term", "B_C_term")
  d_A_N_term_to_A_C_term = cmd.get_distance("A_N_term", "A_C_term")
  d_A_C_term_to_B_N_term = cmd.get_distance("A_C_term", "B_N_term")
  d_A_C_term_to_B_C_term = cmd.get_distance("A_C_term", "B_C_term")
  d_B_N_term_to_B_C_term = cmd.get_distance("B_N_term", "B_C_term")
  #return d_A_N_term_to_B_N_term, d_A_N_term_to_B_C_term, d_A_N_term_to_A_C_term, d_A_C_term_to_B_N_term, d_A_C_term_to_B_C_term, d_B_N_term_to_B_C_term
  to_return = {'d_A_N_term_to_B_N_term' : d_A_N_term_to_B_N_term,
               'd_A_N_term_to_B_C_term' : d_A_N_term_to_B_C_term,
               'd_A_N_term_to_A_C_term' : d_A_N_term_to_A_C_term,
               'd_A_C_term_to_B_N_term' : d_A_C_term_to_B_N_term,
               'd_A_C_term_to_B_C_term' : d_A_C_term_to_B_C_term,
               'd_B_N_term_to_B_C_term' : d_B_N_term_to_B_C_term}
  return to_return



