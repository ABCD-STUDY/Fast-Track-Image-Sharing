#!/bin/bash
#
# Call the anonymizer.sh script for TGZ files containing DICOMs in a directory.
#
# This script will attempt to create a TGZ file from a directory location if
# GE is indicated on the command line (offline multi-band reconstructions).
#
# This file will only work if its run inside the DAIC. Directory locations are
# backed into this script. This script is provided for illustration purposes only.
#


# read in the list of participants from the fast track file
if [ $# -le 2 ]; then
    echo "usage: run_anon.sh fasttrack.json site=(ucsd, chla, ...) scanner=(GE,SI,PH)"
    exit
fi
fasttrack="$1"
site="$2"
option="$3"

if [ ! -d /space/syn09/1/data/mdcornejo/fast-track/$site/ ]; then
    echo "create directory for this site ($site)"
    mkdir /space/syn09/1/data/mdcornejo/fast-track/$site
fi

echo "run anonymizer for site $site, using subjects in $fasttrack"

#get list of participants
participants=`jq -r ".[].pGUID" "$fasttrack" | sort | uniq`

for part in $participants
do
   # tgz files can be in one of two locations
   fl=`ls /space/syn07/1/data/ABCD/incoming/$site/$part*.tgz 2> /dev/null`
   f2=`ls /space/syn08/1/data/ABCD/incoming/$site/$part*.tgz 2> /dev/null`

   # check if we can merge the fl and f2 arrays of file names, maybe we need only one of them?
   if [ ! -z "$fl" ] && [ ! -z "$f2" ]; then
      fl=`echo "$fl $f2"`
   fi
   if [ ! -z "$f2" ] && [ -z "$fl" ]; then
      fl=$f2
   fi

   # now check if the merge array is empty (does the newline in between f1 and f2 work?)
   if [ -z "$fl" ]; then
      continue
   fi
   for file in $fl
   do
       if [ $option == "GE" ]; then
          echo "Do something else for GE data - lookup the unpack data and create a local temporary TGZ for the anonymizer"
          bname=`basename ${file}`
          # we need the study instance UID from this series to get the folder right for input
          studyinstanceuid=`jq -r ".StudyInstanceUID" /space/syn05/1/data/MMILDB/DAL_ABCD_QC/unpack/${site}/${bname%.*}/${bname%.*}.json`

          input="/space/syn05/1/data/MMILDB/DAL_ABCD_QC/unpack/${site}/${bname%.*}/${studyinstanceuid}"

          if [ ! -d "${input}" ]; then
             echo "######  DIRECTORY DOES NOT EXIST ON syn05 #######  $input "  
             input="/space/syn06/1/data/MMILDB/DAL_ABCD_QC/unpack/${site}/${bname%.*}"
             if [ ! -d "${input}" ]; then
                echo "######  DIRECTORY DOES NOT EXIST ON syn06 either #######"
                continue
             fi
          fi
          tout="/tmp/${bname}"
          echo "NOW zip ${input} as ${tout}"
          #create the tgz in /tmp/
           if hash pigz 2>/dev/null; then
          	    tar --dereference -cf - "${input}/" | pigz --fast -p 6 > "${tout}"
           else
                GZIP=-1 tar --dereference -cvzf "${tout}" "${input}/"
           fi
           # now run the anonymizer on this new file
           ./anonymizer.sh -i "${tout}" -d $fasttrack -o /space/syn09/1/data/mdcornejo/fast-track/$site/ -m /space/syn09/1/data/mdcornejo/fast-track/$site/ -p /home/mmilrec14/MetaData/DAL_ABCD_QC/DAL_ABCD_QC_merged_eprime.csv
           # delete the temporary file again to not pollute the system
           /bin/rm -f "/tmp/${bname}"
       else
           ./anonymizer.sh -i $file -d $fasttrack -o /space/syn09/1/data/mdcornejo/fast-track/$site/ -m /space/syn09/1/data/mdcornejo/fast-track/$site/ -p /home/mmilrec14/MetaData/DAL_ABCD_QC/DAL_ABCD_QC_merged_eprime.csv 
       fi
   done
done
