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
# This script looks up and access the files through links in subdirs 'data/incoming' and 'data/unpack', in the current directory.
# New links must be added to these subdirs when a new file system is set up to receive more data.
# Written by Hauke Bartsch ...          Modified by Hauke Bartsch and Octavio Ruiz, 2017oct18-28


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


for part in $participants
do
   # Octavio (2017oct18): ------------------------------------------------------------------------
   # tgz files are in several locations, with more locations added as the project continues.
   # These locations are pointed to by links in one directory:
   fl=`ls data/incoming/*/$site/$part*.tgz 2> /dev/null`
   # ------------------------------------------------------------------------ :Octavio (2017oct18)

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

          # Octavio (2017oct28): ------------------------------------------------------------------------
          #
          # Removed a previous explicit two-directory search for unpack files. Now:
          #
          studyinstanceuid=`jq -r ".StudyInstanceUID"  data/unpack/*/${site}/${bname%.*}/*/*.json`
          seriesinstanceuid=`jq -r ".SeriesInstanceUID" data/unpack/*/${site}/${bname%.*}/*/*.json`
          input=`ls -d data/unpack/*/${site}/${bname%.*}/${studyinstanceuid}`
          if [ ! -d "${input}" ]; then
             echo "UNABLE TO FIND files in UNPACK DIRECTORIES, got \"${input}\""
             continue
          fi
          # ------------------------------------------------------------------------ :Octavio (2017oct28)

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
