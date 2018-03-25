#!/bin/bash
export SUBJECTS_DIR=/Volumes/Users/mbkranz/projects/ACT_Freesurfer_MIke2

cd ${SUBJECTS_DIR}
subList=$(ls -d ACT*_1)

echo $subList

for sub in $subList
do
recon-all -s ${sub} -qcache -target fsaverage -no-isrunning
done

