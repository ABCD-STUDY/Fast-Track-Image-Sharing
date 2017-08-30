#!/bin/bash

# read in the list of participants from the fast track file
if [ $# -le 2 ]; then
    echo "usage: run_anon.sh fasttrack.json site=(ucsd, chla, ...) scanner=(GE,SI,PH)"
    exit
fi
fasttrack="$1"
site="$2"
option="$3"

if [ ! -d /fast-track/$site/ ]; then
    echo "create directory for this site ($site)"
    mkdir /fast-track/$site
fi

echo "run anonymizer for site $site, using subjects in $fasttrack"

#get list of participants

participants=`jq -r ".[].pGUID" "$fasttrack" | sort | uniq`

#participants="NDAR_INVZVD13ZMG"
#participants="NDAR_INVDTHJM3Y9" 
#participants="NDAR_INVG0PD8MBV"

for part in $participants
do

   fl=`ls /space/syn07/1/data/ABCD/incoming/$site/$part*.tgz 2> /dev/null`
   f2=`ls /space/syn08/1/data/ABCD/incoming/$site/$part*.tgz 2> /dev/null`


   # check if we can merge the f1 and f2 arrays of file names, maybe we need only one of them?
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
       echo "WORKING ON FILE:"
       echo ${file}
       if [ $option == "GE" ]; then
          echo "GE scanner"
          echo "Lookup the unpack data and create a local temporary TGZ for the anonymizer"
          bname=`basename ${file}`

          # we need the study instance UID from this series to get the folder right for input
          echo "TRYING syn05"
          studyinstanceuid=`jq -r ".StudyInstanceUID" /space/syn05/1/data/MMILDB/DAL_ABCD_QC/unpack/${site}/${bname%.*}/*/*.json`
          if [ $? != "0" ]; then
              studyinstanceuid="THISISNOTASTUDYINSTANCEUID"
          fi
          seriesinstanceuid=`jq -r ".SeriesInstanceUID" /space/syn05/1/data/MMILDB/DAL_ABCD_QC/unpack/${site}/${bname%.*}/*/*.json`
          input="/space/syn05/1/data/MMILDB/DAL_ABCD_QC/unpack/${site}/${bname%.*}/${studyinstanceuid}"

          if [ ! -d "${input}" ]; then
             echo "DIRECTORY DOES NOT EXIST ON syn05 "
             echo "$input" 
             echo "TRYING syn06"
             studyinstanceuid=`jq -r ".StudyInstanceUID" /space/syn06/1/data/MMILDB/DAL_ABCD_QC/unpack/${site}/${bname%.*}/*/*.json`
             if [ $? != "0" ]; then
                 studyinstanceuid="THISISNOTASTUDYINSTANCEUID"
             fi
             seriesinstanceuid=`jq -r ".SeriesInstanceUID" /space/syn06/1/data/MMILDB/DAL_ABCD_QC/unpack/${site}/${bname%.*}/*/*.json`
             input="/space/syn06/1/data/MMILDB/DAL_ABCD_QC/unpack/${site}/${bname%.*}/${studyinstanceuid}"
             if [ ! -d "${input}" ]; then
                echo "DIRECTORY DOES NOT EXIST ON syn06 either"
                continue
             fi
          fi

          inputJSON=""
          # we have to fix the directory to package, we have to remove the json file inside and we have to put back a better json file (with ClassifyType)
          input2=`ls -d ${input}/* | grep -v json`
          if [ ! -e "${input2}" ]; then
              echo "ERROR: sub-directory with DICOMs does not exist! ${input2}"
              continue
          fi
          inputJSON=`ls ${input}/../*${seriesinstanceuid}*.json`
          inputJSON=`/bin/readlink -e "$inputJSON"`
          if [ ! -e "${inputJSON}" ]; then
              echo "ERROR: Could not find JSON file for this series ${inputJSON}"
              continue
          fi
          echo "INPUT2: $input2, inputJSON: $inputJSON"
          input="$input2"

          #tout="/tmp/${bname}"
          #tout="/space/syn09/1/data/mdcornejo/fast-track/tmp/${bname}"
          tout="/tmp/${bname}"
          echo "NOW ZIPPING ${input}"
          echo "as ${tout}"

          if [ ! -f "${tout}" ]; then
            #create the tgz in /tmp/
            if hash pigz 2>/dev/null; then
          	    tar --dereference -cf - "${input}/" "${inputJSON}" | pigz --fast -p 6 > "${tout}"
            else
                GZIP=-1 tar --dereference -cvzf "${tout}" "${input}/" "${inputJSON}"
            fi
            # now run the anonymizer on this new file
            echo "./anonymizer.sh -i \"${tout}\" -d $fasttrack -o /fast-track/$site/ -m /fast-track/$site/ -p /home/abcddaic/MetaData/DAL_ABCD_QC/DAL_ABCD_QC_combined_eprime.csv"

            ./anonymizer.sh -i "${tout}" -d $fasttrack -o /fast-track/$site/ -m /fast-track/$site/ -p /home/abcddaic/MetaData/DAL_ABCD_QC/DAL_ABCD_QC_incoming_eprime_info.csv
            # delete the temporary file again to not pollute the system
            #/bin/rm -f "/tmp/${bname}"
             /usr/bin/truncate -s 0 "/tmp/${bname}"
          else 
            echo "TAR EXISTS, skip to next entry"
          fi

       else
            echo "Philips/Siemens scanner"
            touchfile=/tmp/`/usr/bin/basename ${file}`
            if [ ! -e "$touchfile" ]; then
               touch "${touchfile}"
               ./anonymizer.sh -i $file -d $fasttrack -o /fast-track/$site/ -m /fast-track/$site/ -p /home/abcddaic/MetaData/DAL_ABCD_QC/DAL_ABCD_QC_incoming_eprime_info.csv
            else
               echo "touch file (${touchfile}) exists already, skip to next entry"
            fi
       fi
   done
done
