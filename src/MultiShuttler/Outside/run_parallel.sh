#!/bin/bash

ARCHITECTURES=(
    "3 3 5 5"
    "4 4 3 3"

    )

# Generate commands
commands=()

for arch in "${ARCHITECTURES[@]}"; do
    # Split the arch string into individual components
    set -- $arch
    commands+=("python3 run_parallel.py $1 $2 $3 $4")  # Pass each argument individually
done

parallel ::: "${commands[@]}"