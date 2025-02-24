#!/bin/bash

# Check if an optimization level is provided
if [ -z "$1" ]; then
    echo "Usage: ./run <optimization_level>"
    echo "Example: ./run 3"
    exit 1
fi

# Set the optimization level based on the input
OPT_LEVEL="-O$1"

# Output executable name
EXECUTABLE="main.out"

g++ -ggdb3 $OPT_LEVEL -std=c++23 -Wall -Wextra -pedantic -o $EXECUTABLE main.cc

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful! Optimization level: $OPT_LEVEL"
    echo "Executable: ./$EXECUTABLE"
    echo "Running $EXECUTABLE..."
    time ./$EXECUTABLE
else
    echo "Compilation failed. Check your code for errors." 
fi

rm -r $EXECUTABLE