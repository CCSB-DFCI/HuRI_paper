#!/bin/bash
output_folder="n_contacts"
mkdir ${output_folder}

pdb_common_path=PATH_TO_EXPERIMENTAL_PDB_FILES
i_done=0
for i in `ls -1 ${pdb_common_path}/ `; do
  PATH_TO_ROSA_PROGRAM ${pdb_common_path}/$i A B -t > ${output_folder}/$i.contacts
  i_done=$((i_done+1))
  if ((i_done%1000==0)); then
    echo "$i_done\n"
  fi
done

