#! /usr/bin/bash --login

# Check if the number of arguments is exactly two
# ###############################################
if [ $# -ne 3 ]; then
  echo "Usage: $0 <dir> <fileNumber> <which-machine>"
  echo "<dir> for which the files are in"
  echo "<fileNumber> = number of file"
  echo "<which-machine> = 'csf3'"
fi

current_dir=$(pwd)

# Load Anaconda Python
if [ "$3" = 'csf3' ]; then
  module load apps/binapps/anaconda3/2023.09  # Python 3.11.5
  source activate qv3d-env
fi

# find the file name
# ##################
#home=$(dirname "$1")
if [[ "$2" =~ ^[0-9]+$ ]]; then
  filename=$(bash $current_dir/bashScripts/getFilename.sh $2)  
else
  echo "Incorrect file number"
  exit 1
fi
fileAddress=$1/h5files/"$filename"

##########
# Get the basename (last part of the address)
basename=$(basename "$1")
# Extract the number part using parameter expansion
density="${basename##*-}"
########################

python $current_dir/pyScripts/plot-contours.py $fileAddress $density

if [ "$3" = 'csf3' ]; then
  conda deactivate
  module unload apps/binapps/anaconda3/2023.09  # Python 3.11.5
fi

#echo "Mission Completed!"
