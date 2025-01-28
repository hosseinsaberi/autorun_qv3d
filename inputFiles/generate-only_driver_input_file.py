#!/usr/bin/env python
"""
Description:
------------
    This script removes the witness beam from the two-beam input file
    to generate an input file only including driver beam.

Created: January 22, 2025
--------

Usage:
------
    python generate-only_driver_input_file.py <input_file(=febe_acc.txt)> <output_file(febe_driver.txt)>")
"""

from datetime import date, datetime
import os
import sys
import shutil


def generat_onlyDriver_inputFile(input_file, output_file):
    # Copy and rename the file
    temp_file = 'temporary.txt'
    shutil.copy(input_file, temp_file)

    with open(temp_file, "r") as file:
        file_contents = file.read()

        # delet lines for & witness beam: Specie2 before &Synchrotron
        # ###########################################################
        # Split content into lines
        lines = file_contents.splitlines()

        # Initialize flags and indices
        start_index = None
        end_index = None

        # Identify the indices for deletion
        for i, line in enumerate(lines):
            if line.startswith("& witness beam: Specie2") and start_index is None:
                start_index = i
            elif line.startswith("&Synchrotron") and start_index is not None:
                end_index = i
                break

        # Rebuild the content by excluding the target lines
        if start_index is not None and end_index is not None:
            lines_to_keep = lines[:start_index] + lines[end_index:]
            file_contents = "\n".join(lines_to_keep)
        else:
            print("ERROR 1- Specified pattern not found in file contents.")

        # delet lines for & witness beam: Specie2 before ######
        # ###########################################################
        # Split content into lines
        lines = file_contents.splitlines()

        # Initialize flags and indices
        start_index = None
        end_index = None

        # Identify the indices for deletion
        for i, line in enumerate(lines):
            if line.startswith("& witness beam: Specie2") and start_index is None:
                start_index = i
            elif line.startswith("# ########################################## #") and start_index is not None:
                end_index = i
                break

        # Rebuild the content by excluding the target lines
        if start_index is not None and end_index is not None:
            lines_to_keep = lines[:start_index] + lines[end_index:]
            file_contents = "\n".join(lines_to_keep)
        else:
            print("ERROR 2-Specified pattern not found in file contents.")

        # delet lines for & witness beam: Specie2 before &Synchrotron
        # ###########################################################
        # Split content into lines
        lines = file_contents.splitlines()

        # Initialize flags and indices
        start_index = None
        end_index = None

        # Identify the indices for deletion
        for i, line in enumerate(lines):
            if line.startswith("& Specie2") and start_index is None:
                start_index = i
            elif line.startswith("&Movie2dHDF5			# save in 2d") and start_index is not None:
                end_index = i
                break

        # Rebuild the content by excluding the target lines
        if start_index is not None and end_index is not None:
            lines_to_keep = lines[:start_index] + lines[end_index:]
            file_contents = "\n".join(lines_to_keep)
        else:
            print("ERROR 3- Specified pattern not found in file contents.")

        # Replace lines
        # #############
        # Loop through the lines and replace the one with "About_inputs"
        lines = file_contents.splitlines()
        for i, line in enumerate(lines):
            if "About_inputs" in line:
                lines[i] = "  About_inputs = Driver simulation (to study energy depletion; wakefield ...) with CLARA FEBE parameters. # short explanation about the script\n"  # Replace the whole line with "aaaaa"

        for i, line in enumerate(lines):
            if "Nspecies" in line:
                lines[i] = "  Nspecies = 2			# includes plasma, driver\n"  # Replace the whole line with "??"

        for i, line in enumerate(lines):
            if "NMovie2dFramesH5" in line:
                lines[i] = "  NMovie2dFramesH5 = 5		# 4, 5, 6 variable outputs in 2d slices\n"  # Replace the whole line with "??"

        for i, line in enumerate(lines):
            if "Frame2" in line:
                lines[i] = " "  # Replace the whole line with "??"
        for i, line in enumerate(lines):
            if "SkipSaveFlag1" in line:
                lines[i] = "  SkipSaveFlag1 = 0		# (=1) to turn off beam particles saving"  # Replace the whole line with "??"
                
        # Join the lines back into a single string
        file_contents = "\n".join(lines)


    with open(output_file, "w") as file:
        file.write(file_contents)

    os.remove(temp_file)





def main(input_file, output_file):
    generat_onlyDriver_inputFile(input_file, output_file)


if __name__ == "__main__":
    # Check if three arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python generate-only_driver_input_file.py <input_file(=febe_acc.txt)> <output_file(febe_driver.txt)>")
    else:
        # Parse command-line arguments
        input_file = sys.argv[1]
        output_file = sys.argv[2]

        # Call the main function with the arguments
        main(input_file, output_file)
