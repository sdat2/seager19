files=( "*.tios" "om_diag.tr" "om_diag.log" "om_diag" "om_month" "om_month.tr" "om_month.log" "om_spin" "om_spin.tr" "om_spin.log" )

old_folder="ocean/SRC/"
new_folder="ocean/old-inputs/"

mkdir ${new_folder}


for i in "${files[@]}"
do
  printf "${old_folder}$i \t ${new_folder}  \n"
  mv ${old_folder}$i ${new_folder}

done
