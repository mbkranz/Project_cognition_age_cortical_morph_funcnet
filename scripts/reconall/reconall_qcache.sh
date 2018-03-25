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

recon-all -s ${sub} -qcache -target fsaverage
echo {$sub}>>$SUBJECTS_DIR/subjectsprocessed_qcache.txt



