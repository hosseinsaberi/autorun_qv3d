#!/bin/bash

# Function to add leading zeros if needed
add_leading_zeros() {
    local num="$1"
    local length=${#num}
    local zeros_needed=$((5 - length))

    if [ $length -lt 5 ]; then
        printf "%0${zeros_needed}d%s\n" 0 "$num"
    else
        echo "$num"
    fi
}

# Check if an argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <integer>"
    exit 1
fi

# Remove leading zeros if present
input_num=$(echo "$1" | sed 's/^0*//')

result=$(add_leading_zeros "$input_num")


# Replace the number in the filename
filename="v2d_mframe_${result}.h5"
echo $filename
