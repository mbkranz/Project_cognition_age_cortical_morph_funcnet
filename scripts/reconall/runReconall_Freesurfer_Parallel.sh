#example parallel call
analdir=/Volumes/HealthAging/ACT/Freesurfer
SUBJECTS_DIR=${analdir}
cd $analdir
allsubs=$(ls -d ACT*_1)

parallel --joblog ${analdir}/logs/log_processRS_$(date +%Y-%m-%d:%H:%M:%S) -j8 --progress --header : ${analdir}/scripts/reconall/reconall_all.sh -s {f1} ::: f1 $(echo $allsubs)

parallel --joblog ${analdir}/logs/log_processRS_$(date +%Y-%m-%d:%H:%M:%S) -j3 --progress --header : ${analdir}/scripts/reconall/reconall_all.sh -s {f1} ::: f1 $(echo $allsubs


#for sub in $allsubs; do
#    recon-all -s ${sub} -qcache -target fsaverage -no-isrunning
#done


parallel --joblog ${analdir}/scripts/logs/log_processRS_$(date +%Y-%m-%d:%H:%M:%S) -j3 --progress --header : ${analdir}/scripts/reconall/reconall_qcache.sh -s {f1} ::: f1 $(echo $allsubs)