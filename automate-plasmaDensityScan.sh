#! /usr/bin/bash

# Check if correct arguments are provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <project-name> <machine-name>"
  echo "project-name = 'FEBE'; 'AWAKE';"
  echo "machinename= 'csf3'  'workstation' 'scarf'"
  exit 1
fi

# Assign command-line arguments to variables
main_directory=$(pwd)

if [ "$1" = 'FEBE-ACC' ]; then
  inputFile="${main_directory}/inputFiles/febe_acc.txt"
  inputName="febe_acc.txt"
elif [ "$1" = 'FEBE-DRIVER' ]; then
  inputFile="${main_directory}/inputFiles/febe_driver.txt"
  inputName="febe_driver.txt"
elif [ "$1" = 'AWAKE' ]; then
  inputFile="${main_directory}/inputFiles/awake_baseline.txt"
  inputName="awake_baseline.txt"
fi

# If the file does not exist, exit
if [ ! -f "$inputFile" ]; then
  echo "File '$inputFile' not found."
  exit 1
fi

# Load Anaconda Python
# ####################
if [ "$2" = "csf3" ]; then
  module load apps/binapps/anaconda3/2023.09  # Python 3.11.5
  source activate qv3d-env
  qv3d_path=""${main_directory}/bin/csf3/qv3dMPIX.e""
elif [ "$2" = "scarf" ]; then
  qv3d_path=""${main_directory}/bin/scarf/qv3dMPIX.e""
else
  qv3d_path=""${main_directory}/bin/workstation/qv3dMPIX.e""
fi

# Assign command-line arguments to variables
if [ "$2" = "scarf" ]; then
  out_directory="/home/vol06/scarf1426/scratch"
else
  out_directory=$(pwd)
fi

# Density scan
# ############
folder_path=${out_directory}/simulation-densScan
# If the folder exists, delete it
if [ -d "${folder_path}" ]; then
  rm -rf "${folder_path}"
fi
mkdir "$folder_path"

word_to_find="plasma_density"
##############################################################################
##############################################################################
# This is the major scan; unblock if needed.
## Define powers array
powers=("e14" "e15" "e16" "e17")
lower=5  # means 5e14
upper=1  # means 1e17
step=0.1
max=9.9

ensity_scan=()

# Loop through the powers and generate values between 8e14 and 5e16 with step 0.1
for power in "${powers[@]}"; do
  if [[ "$power" == "e14" ]]; then
    # Start from 5e14 to 9.9e14
    i=$lower
    while (( $(echo "$i <= $max" | bc -l) )); do
      density_scan+=("$i$power")
      i=$(echo "$i + $step" | bc)
    done
  elif [[ "$power" == "e15" ]]; then
    # Start from 1e16 to 9e16
    i=1
    while (( $(echo "$i <= $max" | bc -l) )); do
      density_scan+=("$i$power")
      i=$(echo "$i + $step" | bc)
    done
  elif [[ "$power" == "e16" ]]; then
    # Start from 1e16 to $upper
    i=1
    while (( $(echo "$i <= $max" | bc -l) )); do
      density_scan+=("$i$power")
      i=$(echo "$i + $step" | bc)
    done
  elif [[ "$power" == "e17" ]]; then
    # Start from 1e16 to $upper
    i=1
    while (( $(echo "$i <= $upper" | bc -l) )); do
      density_scan+=("$i$power")
      i=$(echo "$i + $step" | bc)
    done
  fi
done
##################################################################################
##############################################################################
# Define powers array
#powers=("e15" "e16")
#lower=6  # means 5e14
#upper=3  # means 1e17
#step=0.1
#max=9.9

#ensity_scan=()

# Loop through the powers and generate values between 8e14 and 5e16 with step 0.1
#for power in "${powers[@]}"; do
#  if [[ "$power" == "e15" ]]; then
#    # Start from 5e14 to 9.9e14
#    i=$lower
#    while (( $(echo "$i <= $max" | bc -l) )); do
#      density_scan+=("$i$power")
#      i=$(echo "$i + $step" | bc)
#    done
#  elif [[ "$power" == "e16" ]]; then
#    # Start from 1e16 to 9e16
#    i=1
#    while (( $(echo "$i <= $upper" | bc -l) )); do
#      density_scan+=("$i$power")
#      i=$(echo "$i + $step" | bc)
#    done
#  fi
#done

echo 'plasma density scans over: ' 
echo ' '
echo ${density_scan[@]}
echo 
echo 
##############################################################################
for item in "${density_scan[@]}"; do
  out_dir="$folder_path/density-$item"	# create directory for the new script
  mkdir "$out_dir"
  cp "$inputFile" "${out_dir}"
  cp "${qv3d_path}" "${out_dir}"
  cp "${main_directory}/jobScripts/run_qv3d_on_scarf.sh" "${out_dir}"
  cp "${main_directory}/pyScripts/generate-input-deck.py" "${out_dir}"
  cd "$out_dir"

  # Replace with the new line in the parameters file
  # ************************************************
  new_line="  plasma_density = $item # /cm3"
  sed -i "/\b$word_to_find\b/c\\$new_line" "$inputName"

  python generate-input-deck.py "$inputName"
  rm generate-input-deck.py #"$inputName"

  # Submit job on CSF3
  # ##################
  if [ "$2" = "csf3" ]; then
    # Replace in the job script
    # ************************
    new_line="#$ -N "Density-$item"              # Set the job's name"
    replace='qv3d'
    sed -i "/\b$replace\b/c\\$new_line" jobscript-parallel.sh

    job_id=$(qsub jobscript-parallel.txt)
    cd "$main_directory"

    current_datetime=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$job_id"
    echo "$job_id" at $current_datetime  >> job_IDs.txt
  elif [ "$2" = "scarf" ]; then
    # Replace in the job script
    # ************************
#    new_line="#$ -N "Density-$item"              # Set the job's name"
#    replace='qv3d'
#    sed -i "/\b$replace\b/c\\$new_line" jobscript-parallel.sh

    #job_id=$(sbatch run_qv3d_on_scarf.sh)
    job_id=$(sbatch run_qv3d_on_scarf.sh | awk '{print $4}')
    cd "$main_directory"

    current_datetime=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$job_id"
    #echo "$job_id" at $current_datetime  >> job_IDs.txt	 
    # clear directory
    # ###############
    find . -type f -name '*.dat' -delete
  else
    echo -n "Density $item is running ..."
    mpirun -np 32 ./qv3dMPIX.e v.ini > std.out 2> std.err
    echo 
  fi
  cd ${out_directory}
done

# on CSF3
# #######
if [ "$2" = "csf3" ]; then
  mv job_IDs.txt "$folder_path"
  conda deactivate
  module unload apps/binapps/anaconda3/2023.09  # Python 3.11.5
elif [ "$2" = "scarf" ]; then
  mv job_IDs.txt "$folder_path"
fi

#conda deactivate
#module unload apps/binapps/anaconda3/2023.09  # Python 3.11.5

echo "Mission Completed!"
