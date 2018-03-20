#!/bin/bash
### Usage: 
###     ~]$ sh pattern_Finder.sh "*.c" average_Parameter_Values_nMonths_SHIFT [PATH]
### This will look for average_Parameter_Values_nMonths_SHIFT in files ending with 
### the extension .c
FILES=$1
for f in $3$FILES
do
  echo "Looking for "$2 "in file "$f"..."
  # ls -la $f
  # ls -la $f
  less $f | grep $2
  echo "Please press ENTER to continue. Otherwise Ctrl+C to quit."
  read INPUT_STRING
done
