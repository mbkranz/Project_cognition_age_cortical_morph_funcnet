#!/bin/bash
while getopts “s:p:a:” OPTION
do
  case $OPTION in
    s)
      sub=$OPTARG
      ;;
    ?)
      echo "ERROR: Invalid option"
      #printCommandLine--> from Tim's script (don't know what this does)
      ;;
      esac
 done

recon-all -s ${sub} -all -norm2-b 20 -norm2-n 5 -3T
echo {$sub}>>$SUBJECTS_DIR/subjectsprocessed.txt



