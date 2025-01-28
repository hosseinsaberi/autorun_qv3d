#! /usr/bin/bash

# Check if correct arguments are provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <project-name> <machine-name>"
  echo "project-name = ''AWAKE-BASELINE'; FEBE-ACC'; 'FEBE-DRIVER'"
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

if [ "$2" = "scarf" ]; then
  out_dir="/home/vol06/scarf1426/scratch/simulation"
else
  out_dir="simulation"
fi

# If the folder exists, delete it
for dir in $out_dir*; do
  if [ -d "$dir" ]; then
    rm -rf "${dir}"
  fi
done

# Create run folder - then move files to there
# ############################################
mkdir "$out_dir"
cp "$inputFile" "${out_dir}"
cp "${qv3d_path}" "${out_dir}"
cp "${main_directory}/jobScripts/run_qv3d_on_scarf.sh" "${out_dir}"
cp "${main_directory}/pyScripts/generate-input-deck.py" "${out_dir}"


# in the run folder: generate input deck
# ######################################
cd "$out_dir"
python generate-input-deck.py "$inputName"
rm generate-input-deck.py "$inputName"

# submit the job
# ##############
if [ "$2" = "csf3" ]; then
  qsub jobscript-parallel.txt
  conda deactivate
  module unload apps/binapps/anaconda3/2023.09  # Python 3.11.5
elif [ "$2" = "scarf" ]; then
  sbatch run_qv3d_on_scarf.sh
else
  echo "Running ..."
  mpirun -np 32 ./qv3dMPIX.e v.ini > std.out 2> std.err
fi

# clear directory
# ###############
find . -type f -name '*.dat' -delete

echo "Mission Completed!"
