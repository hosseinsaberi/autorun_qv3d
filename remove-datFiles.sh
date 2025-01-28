#!/bin/bash

# Function to remove *.dat files from current directory and subdirectories
remove_dat_files() {
    # Remove *.dat files in current directory
    find . -type f -name '*.dat' -delete

    # Remove *.dat files in subdirectories
    find . -type d -exec sh -c 'cd "$0" && find . -type f -name "*.dat" -delete' {} \;
}

# Call the function to remove *.dat files
remove_dat_files

