#!/usr/bin/env python3
"""
Anonymize a list of DICOM files

List of tags for example from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3229700/

(This script runs on tgz files to produce a new tgz file with values)

If a *.hdr file is detected inside the input tgz no output is generated.

This script and associated resources are in: https://github.com/ABCD-STUDY/Fast-Track-Image-Sharing

Written by Hauke, Daniela.  Modified by Octavio Ruiz, 2017oct04
"""

import sys, os, time, atexit, stat, tempfile, copy, tarfile, datetime, io, getopt
import dicom, json, re, logging, logging.handlers, threading, pycurl, csv
import struct
from signal import SIGTERM
from dicom.filereader import InvalidDicomError
from dicom.dataelem import DataElement
from dicom.tag import Tag
from dicom.filebase import DicomBytesIO
from multiprocessing import Pool, TimeoutError
from io import StringIO
import sqlite3
import pandas as pd
import hashlib

"""
  Do the main work here
"""
def anonymize(df, anonInfo, mode="dir", tarout=None):
        try:
                # read the file using pydicom
                dataset = dicom.read_file(df)
        except IOError:
                print("Error: reading dicom file %s" % df)
                log.error("Error: reading dicom file %s" % df)
                return 0
        except OSError:
                print("Error: reading dicom file %s" % df)
                log.error("Error: reading dicom file %s" % df)
                return 0
        except InvalidDicomError:
                print("Warning: not a dicom file %s" % df)
                log.info("Warning: not a dicom file %s" % df)
                # maybe this is not a dicom file but the json that contains series information,
                # lets add the file to the output as well.
                with open(df) as data_file:
                    try:
                        info = json.load(data_file)
                    except:
                        print("Error: found file neither DICOM nor JSON (%s, %s) in %s" % (anonInfo['origFileName'], sys.exc_info()[0], anonInfo['inputfile']))
                        log.error("Error: found file neither DICOM nor JSON")
                        pass

                    try:
                        # replace information from the json file
                        info['PatientName'] =  ['pGUID']
                        info['PatientID'] = anonInfo['pGUID']
                        info['StudyDescription'] = "Adolescent Brain Cognitive Development Study"
                        info['SeriesDescription'] = ', '.join(anonInfo['ClassifyType'])
                        now = datetime.datetime.now()
                        info['Anonymized'] = now.strftime("fast track v1.0, @DAIC %B %Y")
                        info['NumFiles'] = "removed"
                        info['IncomingConnection'] = "removed"
                        bday = datetime.datetime.strptime(anonInfo['dob'], '%Y-%m-%d')
                        info['DateOfBirth'] = bday.strftime("%Y%m%d")
                        anonInfo['DateOfBirth'] = info['DateOfBirth']
                        sday = datetime.datetime.strptime(anonInfo['StudyDate'], '%Y%m%d')
                        info['PatientsAge'] = ("%.1f" % ((sday-bday).days/365.25))   # PatientAge
                        anonInfo['PatientsAge'] = info['PatientsAge']
                        info['PatientsSex'] = anonInfo['gender']
                        anonInfo['RepetitionTime'] = '0'
                        try: 
                                anonInfo['RepetitionTime'] = info['RepetitionTime']
                        except:
                                pass
                        anonInfo['EchoTime'] = '0'
                        try:
                                anonInfo['EchoTime'] = info['EchoTime']
                        except:
                                pass
                        # calculate a series time without the fractional seconds
                        seriestime = anonInfo['SeriesTime']
                        if "." in seriestime:
                                seriestime = seriestime.split(".")[0]
                        bidstype = bidsFolderByClassifyType( anonInfo['ClassifyType'] )

                        if 'ABCD-SST-fMRI' in anonInfo['ClassifyType'] or 'ABCD-MID-fMRI' in anonInfo['ClassifyType'] or 'ABCD-nBack-fMRI' in anonInfo['ClassifyType']: 
                                # lookup the EPrime file for SeriesInstanceUID    
                                eventEdatFileName = "sub-%s/ses-%s/%s/%s_run-%s%s-EventRelatedInformation.%s" % (anonInfo['pGUID_BIDS'],
                                                                               anonInfo['event_BIDS'],
                                                                               bidstype,
                                                                               anonInfo['ABCDType'],
                                                                               anonInfo['SeriesDate'],
                                                                               seriestime,
                                                                               eprime_flags[1])
                                if not os.path.isfile(eprime_flags[0]):
                                        print("Error: eprime file not found, path exists in spreadsheet but file is not there\n")
                                else:
                                        tarout.add(eprime_flags[0],eventEdatFileName)

                        # if we have type ABCD-DTI we should add bval and bvec files now
                        if 'ABCD-DTI' in anonInfo['ClassifyType']:
                                subjectName = os.path.basename(tarout.name)
                                subjectName, fext = os.path.splitext(subjectName)

                                GEbvals  = """0 0 3000 3000 2000 3000 1000 3000 3000 2000 3000 1000 3000 3000 500 3000 2000 3000 3000 1000 3000 500 3000 3000 2000 3000 1000 3000 3000 2000 3000 1000 3000 3000 0 3000 2000 3000 1000 3000 3000 500 3000 2000 3000 3000 1000 3000 0 3000 3000 2000 3000 1000 3000 3000 2000 3000 1000 3000 3000 500 3000 2000 3000 3000 1000 3000 0 3000 2000 3000 3000 1000 3000 500 3000 3000 2000 3000 1000 3000 3000 2000 3000 1000 3000 3000 0 3000 2000 3000 3000 1000 3000 500 3000 3000 2000 3000 1000 3000 0"""
                                GEbvects = """0.000000 0.000000 0.654875 0.271924 -0.957387 -0.118810 -0.957387 0.326879 -0.394817 -0.190618 0.906253 -0.190617 0.566330 -0.558540 -0.006227 0.648904 -0.325492 0.435216 -0.060085 0.325492 0.290295 -0.636679 0.008761 0.098843 0.022796 -0.162887 -0.022796 0.095045 -0.214960 -0.596978 0.741627 -0.596977 0.847852 0.943158 -0.000000 0.215650 0.416251 0.060360 0.416250 0.814723 -0.265965 -0.399916 -0.441982 -0.734448 0.566124 -0.716591 0.734448 -0.331272 -0.000000 -0.465735 -0.480716 0.581881 -0.961363 0.581881 -0.920360 0.393331 -0.928822 0.852276 0.928822 -0.571183 -0.535413 0.789559 0.442791 -0.861660 0.143848 0.893460 -0.861660 0.795825 -0.000000 0.399108 0.525770 0.806535 0.691380 -0.525769 0.973847 -0.887690 0.052107 0.696694 -0.034273 0.626272 -0.034274 0.994108 0.738939 -0.084507 0.702182 -0.084509 -0.505053 -0.371230 -0.000000 0.050870 -0.630413 0.170569 0.312160 -0.630413 -0.887018 0.152552 0.605767 -0.215791 -0.491239 0.210754 0.491239 0.787853 -0.000000
0.000000 0.000000 0.355659 0.933965 0.187071 0.110437 0.187070 0.866547 0.108337 0.585245 0.349334 0.585244 0.470917 -0.665203 0.064446 -0.611989 -0.844583 0.497941 -0.972644 0.844583 0.619397 -0.653135 -0.548745 -0.637960 0.678409 0.917294 -0.678409 0.821836 -0.721063 0.258723 -0.039918 0.258722 -0.200363 0.304222 0.000000 -0.732737 0.343103 -0.223624 0.343104 -0.561120 -0.335097 0.828420 -0.764517 0.292307 0.807614 -0.584614 -0.292308 -0.943028 0.000000 0.168038 0.789299 0.362502 0.253123 0.362503 -0.003999 -0.915435 0.354407 -0.324551 -0.354407 0.799445 -0.167226 -0.384930 0.203222 -0.322776 0.431269 -0.414910 -0.322776 0.511407 0.000000 -0.523871 -0.761172 0.590581 -0.710203 0.761173 -0.080055 -0.101313 -0.996034 0.034360 -0.991162 -0.286200 -0.991162 0.072683 -0.519516 -0.056169 0.353175 -0.056169 0.581614 0.770180 0.000000 -0.858370 -0.756662 -0.325484 -0.422488 -0.756662 -0.152189 0.851204 0.776332 -0.937112 0.826626 0.013433 -0.826626 0.274495 0.000000
0.000000 0.000000 0.666817 0.231877 -0.220034 0.986756 -0.220034 0.377155 -0.912350 0.788133 -0.238058 0.788133 -0.676393 -0.495518 -0.997902 -0.452098 -0.425129 0.750095 -0.224396 0.425129 -0.729435 0.409944 0.835944 -0.763700 -0.734331 -0.363371 0.734330 -0.561740 -0.658681 -0.759395 -0.669624 -0.759395 0.490920 0.133799 0.000000 0.645439 0.842031 0.972805 0.842031 -0.146185 0.903865 -0.392156 0.469218 0.612489 -0.165116 0.380425 -0.612489 -0.030940 0.000000 0.868823 -0.381994 -0.728016 0.108216 -0.728016 0.391053 0.085260 -0.108098 -0.410233 0.108097 0.186110 0.827870 -0.477939 0.873291 -0.391608 0.890682 0.171983 -0.391608 0.324230 0.000000 -0.752511 0.379714 -0.026762 0.132686 -0.379714 0.212633 -0.449158 0.072118 0.716545 0.128152 0.725171 0.128153 -0.080414 0.429037 0.994838 0.618230 0.994838 -0.637689 0.518664 0.000000 -0.510502 0.173326 -0.930036 0.850917 0.173326 -0.435932 0.502174 0.174228 0.274326 0.274543 -0.977447 -0.274542 -0.551307 0.000000"""
                                Philipsbvals  = """0 1000 1000 1000 1000 1000 1000 1000 1000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 2000 2000 2000 2000 2000 2000 2000 0 0 0 500 500"""
                                Philipsbvects = """-0.000000 0.597038 0.861645 -0.957395 -0.734424 0.525812 0.190702 -0.416245 -0.491185 0.265989 0.371191 -0.691391 -0.852281 0.654887 -0.143797 -0.973848 0.215804 -0.170604 -0.696719 -0.741627 -0.399098 0.480702 -0.326906 -0.795827 -0.095005 -0.566311 0.920349 -0.806513 -0.393316 -0.060400 -0.050899 0.465720 -0.943169 -0.008801 -0.566116 0.571216 0.394822 -0.994108 -0.052102 -0.581885 0.034292 0.325535 0.861645 0.190567 -0.416245 0.630401 -0.957395 -0.943169 -0.491185 -0.789577 0.399988
-0.000000 -0.258769 0.322835 -0.187079 0.292384 -0.761181 -0.585269 -0.343146 0.826667 0.335086 -0.770182 0.710190 0.324593 -0.355693 -0.431290 0.080104 0.937116 0.325507 -0.034401 0.039901 0.523898 -0.789304 -0.866515 -0.511417 -0.821846 -0.470910 0.004000 -0.590609 0.915437 0.223602 0.858378 -0.168007 -0.304190 0.548733 -0.807622 -0.799423 -0.108306 -0.072701 0.996035 -0.362530 0.991167 0.844578 0.322835 -0.585294 -0.343146 0.756677 -0.187079 -0.304190 0.826667 0.384867 -0.828389
0.000000 -0.759331 -0.391594 -0.219991 -0.612482 -0.379638 0.788095 0.842016 -0.274516 0.903862 0.518688 0.132698 -0.410191 0.666787 0.890680 0.212610 0.274305 -0.930021 0.716519 -0.669625 -0.752497 -0.382002 0.377206 0.324211 -0.561732 -0.676414 0.391078 -0.026800 0.085303 0.972807 -0.510487 0.868837 0.133796 0.835951 -0.165105 0.186105 -0.912351 -0.080401 0.072103 -0.727999 0.128106 -0.425106 -0.391594 0.788108 0.842016 0.173308 -0.219991 0.133796 -0.274516 -0.477960 -0.392150"""
                                Philipsbvals2 = """0 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 500 500 500 500 1000 1000 1000 1000 1000 1000 1000 0 0 2000 2000 2000 2000 2000 2000 2000 2000"""
                                Philipsbvects2 = """-0.000000 -0.312199 -0.215709 -0.435210 0.961373 -0.605806 0.558519 -0.893470 0.214990 0.505088 0.442006 -0.098798 0.535393 0.331306 -0.847867 0.271891 -0.648897 -0.210808 -0.442797 -0.738939 0.887025 -0.787873 0.716608 -0.702198 -0.906263 -0.814723 0.162897 -0.118795 0.060102 -0.290307 -0.152605 0.887651 0.636661 0.006124 -0.928844 0.084520 -0.581836 0.034297 0.022862 0.630406 -0.325478 0.214990 -0.702198 0.084520 0.596987 -0.957385 0.734451 -0.022780 0.491256 0.928809 -0.525785
-0.000000 0.422499 0.732730 -0.497912 -0.253093 -0.776307 0.665222 0.414886 0.721065 -0.581587 0.764510 0.637986 0.167198 0.943017 0.200392 -0.933970 0.611997 -0.013400 -0.203198 0.519527 0.152204 -0.274491 0.584607 -0.353199 -0.349286 0.561116 -0.917284 -0.110396 0.972640 -0.619414 -0.851208 0.101404 0.653074 -0.064420 0.354357 0.056116 -0.362543 0.991157 0.678398 0.756660 -0.844614 0.721065 -0.353199 0.056116 -0.258686 -0.187019 -0.292335 -0.678384 -0.826599 -0.354427 0.761182
0.000000 0.850897 0.645427 0.750117 0.108197 0.174202 -0.495516 0.171994 -0.658668 -0.637685 0.469206 -0.763684 0.827889 -0.030901 0.490881 0.231893 -0.452098 -0.977436 0.873293 0.429023 -0.435912 -0.551281 0.380404 0.618198 -0.238090 -0.146204 -0.363394 0.986763 -0.224409 -0.729416 0.502152 -0.449214 0.410070 -0.997904 0.108074 0.994840 -0.728031 0.128182 0.734339 0.173362 0.425078 -0.658668 0.618198 0.994840 -0.759400 -0.220087 0.612471 -0.734354 0.274594 -0.108141 0.379673"""
                                Siemensbvals = """0 0 3000 3000 2000 3000 1000 3000 3000 2000 3000 1000 3000 3000 500 3000 2000 3000 3000 1000 3000 500 3000 3000 2000 3000 1000 3000 3000 2000 3000 1000 3000 3000 0 3000 2000 3000 1000 3000 3000 500 3000 2000 3000 3000 1000 3000 0 3000 3000 2000 3000 1000 3000 3000 2000 3000 1000 3000 3000 500 3000 2000 3000 3000 1000 3000 0 3000 2000 3000 3000 1000 3000 500 3000 3000 2000 3000 1000 3000 3000 2000 3000 1000 3000 3000 0 3000 2000 3000 3000 1000 3000 500 3000 3000 2000 3000 1000 3000 0"""
                                Siemensbvects = """-0.000000 -0.000000 0.355858 0.934239 0.186957 0.110723 0.186895 -0.866933 -0.108187 -0.585872 -0.349306 -0.586346 -0.470522 0.664792 -0.063733 0.611729 0.844008 -0.498377 0.972395 -0.845213 -0.618837 0.653223 0.549150 0.637402 -0.677539 -0.916902 0.679434 -0.821220 0.720519 -0.258598 0.039868 -0.258004 0.200443 -0.304328 0.000000 0.733334 -0.343553 0.223929 -0.343535 0.561079 0.335306 -0.827192 0.764866 -0.292326 -0.807523 0.584813 0.292317 0.943019 0.000000 -0.168135 -0.788922 -0.362235 -0.253089 -0.361941 0.004031 0.915551 -0.354290 0.324479 0.354392 -0.799542 0.167535 0.384375 -0.203540 0.322526 -0.431705 0.414904 0.322511 -0.511701 0.000000 0.523460 0.761599 -0.590579 0.710286 -0.760425 0.080122 0.101249 0.996112 -0.034504 0.991356 0.286578 0.991413 -0.072758 0.519933 0.056411 -0.353456 0.056509 -0.581172 -0.770686 0.000000 0.857754 0.756794 0.325031 0.422880 0.757200 0.152176 -0.852380 -0.776488 0.937393 -0.827044 -0.013243 0.826057 -0.274353 0.000000
-0.000000 -0.000000 -0.655513 -0.272029 0.957094 0.118828 0.956911 0.327157 -0.394373 -0.190835 0.905992 -0.191189 0.565902 -0.558312 -0.005951 0.648472 -0.325327 0.435523 -0.059997 0.326103 0.289867 -0.637838 0.008708 0.098754 0.022966 -0.162769 -0.022751 0.094978 -0.214940 -0.596246 0.741044 -0.595893 0.848338 0.943300 -0.000000 0.215886 0.416952 0.060414 0.417023 0.814552 -0.266120 -0.399771 -0.442400 -0.735112 0.565939 -0.716992 0.733233 -0.331180 -0.000000 -0.466156 -0.480470 0.581260 -0.961493 0.581147 -0.920818 0.393314 -0.928715 0.851797 0.928998 -0.571365 -0.535979 0.788205 0.443093 -0.861221 0.144001 0.893658 -0.860948 0.796077 -0.000000 0.398642 0.525984 0.806489 0.691512 -0.525146 0.974071 -0.886396 0.052210 0.697290 -0.034320 0.626752 -0.034278 0.994016 0.739273 -0.084637 0.702679 -0.084749 -0.504541 -0.371331 -0.000000 0.050827 -0.630628 0.170539 0.312558 -0.630249 -0.886549 0.153086 0.605923 -0.215958 -0.491445 0.210576 0.490836 0.787291 -0.000000
0.000000 0.000000 0.666084 0.230647 -0.221399 0.986722 -0.222244 0.376026 -0.912560 0.787614 -0.239087 0.787175 -0.677026 -0.496326 -0.997949 -0.453070 -0.426396 0.749627 -0.225495 0.423405 -0.730081 0.407998 0.835678 -0.764177 -0.735128 -0.364413 0.733384 -0.562651 -0.659282 -0.760011 -0.670272 -0.760490 0.490046 0.132552 0.000000 0.644682 0.841500 0.972731 0.841472 -0.147297 0.903742 -0.394888 0.468254 0.611683 -0.166191 0.379362 -0.613938 -0.032168 0.000000 0.868579 -0.383080 -0.728644 0.107132 -0.728880 0.389973 0.084088 -0.109398 -0.411284 0.106625 0.185135 0.827441 -0.480613 0.873063 -0.392779 0.890446 0.170968 -0.393388 0.323145 0.000000 -0.753044 0.378560 -0.028144 0.131548 -0.382067 0.211578 -0.451719 0.070955 0.715958 0.126633 0.724607 0.126198 -0.081473 0.427954 0.994814 0.617504 0.994799 -0.638496 0.517838 0.000000 -0.511541 0.171962 -0.930200 0.850576 0.171567 -0.436891 0.500012 0.172985 0.273234 0.272912 -0.977488 -0.276964 -0.552180 0.000000"""

                                bvalName = "sub-%s/ses-%s/%s/%s_run-%s%s.bval" % (anonInfo['pGUID_BIDS'],
                                                                               anonInfo['event_BIDS'],
                                                                               bidstype,
                                                                               anonInfo['ABCDType'],
                                                                               anonInfo['SeriesDate'],
                                                                               seriestime)
                                bvecName = "sub-%s/ses-%s/%s/%s_run-%s%s.bvec" % (anonInfo['pGUID_BIDS'],
                                                                               anonInfo['event_BIDS'],
                                                                               bidstype,
                                                                               anonInfo['ABCDType'],
                                                                               anonInfo['SeriesDate'],
                                                                               seriestime)
                                if 'GE' in anonInfo['ClassifyType']:
                                        tinfo = tarfile.TarInfo(name=bvalName)
                                        tinfo.size = len(GEbvals)
                                        tarout.addfile(tinfo, io.BytesIO(GEbvals.encode('utf8')))
                                        tinfo = tarfile.TarInfo(name=bvecName)
                                        tinfo.size = len(GEbvects)
                                        tarout.addfile(tinfo, io.BytesIO(GEbvects.encode('utf8')))
                                elif 'PHILIPS' in anonInfo['ClassifyType']:
                                        m  = re.compile('DTI\s*([0-9]+)')
                                        mm = m.match(anonInfo['SeriesDescription'])
                                        #log.error("now run ...%s \n" % mm.group(1)) 
                                        if mm and mm.group(1) == '1':
                                                #log.error("Detected Philips DTI 1")
                                                tinfo = tarfile.TarInfo(name=bvalName)
                                                tinfo.size = len(Philipsbvals)
                                                tarout.addfile(tinfo, io.BytesIO(Philipsbvals.encode('utf8')))
                                                tinfo = tarfile.TarInfo(name=bvecName)
                                                tinfo.size = len(Philipsbvects)
                                                tarout.addfile(tinfo, io.BytesIO(Philipsbvects.encode('utf8')))
                                        else:
                                                #log.error("Detected Philips DTI 2")
                                                bvalName2 = "%s" % bvalName
                                                bvecName2 = "%s" % bvecName
                                                tinfo = tarfile.TarInfo(name=bvalName2)
                                                tinfo.size = len(Philipsbvals2)
                                                tarout.addfile(tinfo, io.BytesIO(Philipsbvals2.encode('utf8')))
                                                tinfo = tarfile.TarInfo(name=bvecName2)
                                                tinfo.size = len(Philipsbvects2)
                                                tarout.addfile(tinfo, io.BytesIO(Philipsbvects2.encode('utf8')))
                                if 'SIEMENS' in anonInfo['ClassifyType']:
                                        tinfo = tarfile.TarInfo(name=bvalName)
                                        tinfo.size = len(Siemensbvals)
                                        tarout.addfile(tinfo, io.BytesIO(Siemensbvals.encode('utf8')))
                                        tinfo = tarfile.TarInfo(name=bvecName)
                                        tinfo.size = len(Siemensbvects)
                                        tarout.addfile(tinfo, io.BytesIO(Siemensbvects.encode('utf8')))                                        
                                        

                        with tempfile.NamedTemporaryFile() as temp:
                                with open(temp.name, 'w') as outfile:
                                    json.dump(info, outfile, indent=4, sort_keys=True)
                                print("JSON encoding detected...\n")
                                log.info("JSON encoding detected...")
                                subjectName = os.path.basename(tarout.name)
                                subjectName, fext = os.path.splitext(subjectName)

                                imageName = "sub-%s/ses-%s/%s/%s_run-%s%s.json" % (anonInfo['pGUID_BIDS'], 
                                                                                   anonInfo['event_BIDS'],
                                                                                   bidstype, 
                                                                                   anonInfo['ABCDType'], 
                                                                                   anonInfo['SeriesDate'], 
                                                                                   seriestime)
                                tinfo = tarfile.TarInfo(name=imageName)
                                tarout.add(temp.name, imageName)
                    except:
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                        print(exc_type, fname, exc_tb.tb_lineno)
                        print("Error: cannot get value of key (%s, %s)" % (anonInfo['origFileName'], sys.exc_info()[0]))
                        log.error("Error: found file neither DICOM nor JSON")
                        pass
                return 0
        num = 0
        for tagEntry in anonymizeThis:
                # remove this tag from this DICOM file
                newValue = ""
                if 'value' in tagEntry:
                        newValue = tagEntry['value']
                if 'group' in tagEntry:
                        if not 'element' in tagEntry:
                                print("Error: anonymizeTag element contains group but no tag")
                                log.error("Error: anonymizeTag element contains group but no tag")
                        else:
                                try:
                                        dataset[int(tagEntry['group'],0),int(tagEntry['element'],0)].value = newValue
                                        num = num + 1
                                except KeyError:
                                        # if we don't find this key in the data group, maybe its from the meta data group
                                        try: 
                                                dataset.file_meta[int(tagEntry['group'],0),int(tagEntry['element'],0)].value = newValue
                                                num = num + 1
                                        except KeyError:
                                                #print "Erro: Key %s:%s does not exist" % (tagEntry['group'], tagEntry['element'])
                                                #try:
                                                #        log.error("Error: Key %s:%s does not exist (meta or data)" % (tagEntry['group'], tagEntry['element']))
                                                #except FileNotFoundError:
                                                #        pass
                                                pass
                                        pass
                if 'name' in tagEntry:
                        #newData = DataElement(dataset.data_element(tagEntry['name']).tag, dataset.data_element(tagEntry['name']).VR, newValue)
                        #dataset.add(newData)
                        try:
                                dataset.data_element(tagEntry['name']).value = newValue
                                num = num + 1
                        except KeyError:
                                try:
                                        dataset.file_meta.data_element(tagEntry['name']).value = newValue
                                        num = num + 1
                                except KeyError:
                                      print("Error: got KeyError trying to set the value of %s" % tagEntry['name'])
                                      #log.error("Error: got KeyError trying to set the value of %s" % tagEntry['name'])
                                      pass
                                pass
        
        dataset[int("0x10",0),int("0x10",0)].value = anonInfo['pGUID']
        dataset[int("0x10",0),int("0x20",0)].value = anonInfo['pGUID']
        dataset[int("0x10",0),int("0x40",0)].value = anonInfo['gender']
        
        
        #hash for DeviceSerialNumber 
        encoded = dataset[int("0x18",0),int("0x1000",0)].value.encode('utf-8')
        h = "anon%s" % hashlib.sha224(encoded).hexdigest()
        dataset[int("0x18",0),int("0x1000",0)].value = h[0:8]
        print("encoded: %s" % encoded)
        print("h: %s" % h[0:8])

        #return this value from dicom if not there already
        if not 'SoftwareVersion' in anonInfo:
            anonInfo['SoftwareVersion'] = ''
            try:
                anonInfo['SoftwareVersion'] = dataset[int("0x18",0),int("0x1020",0)].value
            except KeyError:
                pass
        if not 'FlipAngle' in anonInfo:
            anonInfo['FlipAngle'] = ''
            try:
                anonInfo['FlipAngle'] = dataset[int("0x18",0),int("0x1314",0)].value
            except KeyError:
                pass
        if not 'AcquisitionMatrix' in anonInfo:
            anonInfo['AcquisitionMatrix'] = ''
            try:
                anonInfo['AcquisitionMatrix'] = dataset[int("0x18",0),int("0x1310",0)].value
            except KeyError:
                pass
        if not 'FoV' in anonInfo:
            anonInfo['FoV'] = ''
            try:
                anonInfo['FoV'] = dataset[int("0x51",0),int("0x100c",0)].value
            except KeyError:
                pass
        if not 'PatientPosition' in anonInfo:
            anonInfo['PatientPosition'] = ''
            try:
                anonInfo['PatientPosition'] = dataset[int("0x18",0),int("0x5100",0)].value
            except KeyError:
                pass
        if not 'PhotometricInterpretation' in anonInfo:
            anonInfo['PhotometricInterpretation'] = ''
            try:
                anonInfo['PhotometricInterpretation'] = dataset[int("0x28",0),int("0x04",0)].value
            except KeyError:
                pass
        if not 'ReceiveCoilName' in anonInfo:
            anonInfo['ReceiveCoilName'] = ''
            try:
                anonInfo['ReceiveCoilName'] = dataset[int("0x18",0),int("0x1250",0)].value
            except KeyError:
                pass
        if not 'TransmitCoilName' in anonInfo:
            anonInfo['TransmitCoilName'] = ''
            try:
                anonInfo['TransmitCoilName'] = dataset[int("0x18",0),int("0x1251",0)].value
            except KeyError:
                pass
        if not 'SliceThickness' in anonInfo:
            anonInfo['SliceThickness'] = ''
            try:
                anonInfo['SliceThickness'] = dataset[int("0x18",0),int("0x50",0)].value
            except KeyError:
                pass
        if not 'ImagePixelSpacing' in anonInfo:
            anonInfo['ImagePixelSpacing'] = ''
            try:
                anonInfo['ImagePixelSpacing'] = dataset[int("0x18",0),int("0x1050",0)].value
            except KeyError:
                pass
        # calculate the Age by using dob and scan date
        bday = datetime.datetime.strptime(anonInfo['dob'], '%Y-%m-%d')
        dataset[int("0x10",0),int("0x30",0)].value = bday.strftime('%Y%m%d')
        sday = datetime.datetime.strptime(dataset[int("0x08",0),int("0x20",0)].value, '%Y%m%d')
        dataset[int("0x10",0),int("0x1010",0)].value = ("%.1f" % ((sday-bday).days/365.25))   # PatientAge
        try:
                dataset[int("0x08",0),int("0x1030",0)].value = "Adolescent Brain Cognitive Development Study"
        except KeyError:
                pass
        dataset[int("0x08",0),int("0x103e",0)].value = ''.join([', '.join(sorted(anonInfo['ClassifyType'])), ' (', anonInfo['event'], ')'])
        

        # find out if the file is a symbolic link, overwrite the origin instead
        #if os.path.islink(df):
        #        print("File is symlink, replace origin instead")
        #        log.error("File is symlink %s" % df)
        # and overwrite the file again

        try:
                #print("Save the file again in: %s" % df)
                if mode == "dir":
                    log.error("Save the file again in: %s" % df)
                    dataset.save_as(df)
                else:
                    # calculate a series time without the fractional seconds
                    seriestime = anonInfo['SeriesTime']
                    if "." in seriestime:
                            seriestime = seriestime.split(".")[0]

                    bidstype = bidsFolderByClassifyType( anonInfo['ClassifyType'] )

                    # assumption is that we need to save into a tgz now
                    # this is potentially a critical section, we should only try to write to the tar file if no one else does
                    with tempfile.NamedTemporaryFile() as temp:
                        dataset.save_as(temp.name)
                        subjectName = os.path.basename(tarout.name)
                        subjectName, fext = os.path.splitext(subjectName)
                        imageName = "sub-%s/ses-%s/%s/%s_run-%s%s/sub-%s_ses-%s_dicom%06d.dcm" % (anonInfo['pGUID_BIDS'], anonInfo['event_BIDS'], bidstype,
                                                                                        anonInfo['ABCDType'], anonInfo['SeriesDate'], seriestime, anonInfo['pGUID_BIDS'], 
                                                                                        anonInfo['event_BIDS'], dataset[int("0x20",0),int("0x13",0)].value)
                        tinfo = tarfile.TarInfo(name=imageName)
                        tarout.add(temp.name, imageName)
        except:
                print("Error: could not overwrite DICOM %s %s" % (df, sys.exc_info()[0]))
                log.error("Error: could not overwrite DICOM %s" % df)
        log.info("Anonymize %s (%d tags touched)" % (df, num))
        return 0

def createMetaDataDB( metadatadir, metadata ):
    sqlite_file = ''.join([metadatadir, '/', 'metadata.sqlite'])    # name of the sqlite database file
    table_name1 = 'image03'  # name of the table to be created
    table_name2 = 'iqc'  # name of the table to be created

    # Connecting to the database file
    conn = sqlite3.connect(sqlite_file)
    c = conn.cursor()

    # Creating a new SQLite table
    c.execute('CREATE TABLE {tn} (id INTEGER PRIMARY KEY)'.format(tn=table_name1))
    for key in metadata:
        c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}".format(tn=table_name1,cn=key,ct="TEXT"))

    # Creating a second table with 1 column and set it as PRIMARY KEY
    # note that PRIMARY KEY column must consist of unique values!
    #iqctable = { 'someColumn': 'SomeValue' }
    #c.execute('CREATE TABLE {tn}'.format(tn=table_name2))
    #for key in iqctable:
    #    c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}".format(tn=table_name2,cn=key,ct="TEXT"))

    # Committing changes and closing the connection to the database file
    conn.commit()
    conn.close()


# ---------------------------------------------------------------------------------------------------------
def uploadToNDA( metadatadir, metadata ):

    import requests
    import subprocess

    imagefilename = metadata['image_file']

    # --------------------------------------------------------
    #               Upload metadata to miNDAR

    try:
        with open('login_credentials.json','r') as f:
            try:
                login_credentials = json.load(f)
            except ValueError:
                print("Error: could not read login_credentials.json in the current directory or syntax error")
                log.error("Error: could not read login_credentials.json in the current directory or syntax error")
                sys.exit(0)
    except IOError:
        print("Error: unable to read login_credentials.json file in the current directory")
        log.error("Error: could not read login_credentials.json in the current directory")
        sys.exit(0)
        
    username = login_credentials['miNDAR']['username']
    password = login_credentials['miNDAR']['password']

    # NDA requires a number in the first image_resolution field
    metadata['image_resolution1'] = '0'

    # Assembly package from metadata
    package = {
        "schemaName": "abcd_upload_107927",
        "dataStructureRows": [ {
                "shortName":  "image03",
                "dataElement": []
        } ]
    }
    for i,v in metadata.items():
        t = v
        if isinstance(t, list):
            t = json.dumps(t)
        package['dataStructureRows'][0]['dataElement'].append( { "name": i, "value": t } )


    print( json.dumps(package, indent=2) )
    input('ENTER to continue')

    # Upload metadata package
    res = requests.post( "https://ndar.nih.gov/api/mindar/import",
                        auth=requests.auth.HTTPBasicAuth(username, password),
                        headers={'content-type':'application/json'},
                        data = json.dumps(package) )
    miNDA_ok  = res.ok
    miNDA_msg = res.text
    # --------------------------------------------------------

    # --------------------------------------------------------
    #                  Upload file to aws

    rs  =  subprocess.run( ['/home/oruiz/.local/bin/aws', 's3', 'cp', imagefilename, 's3://nda-abcd/'], stderr=subprocess.PIPE )
    S3_ok  =  (rs.returncode == 0)
    S3_msg =  rs.stderr
    # --------------------------------------------------------
    
    return [miNDA_ok, miNDA_msg, S3_ok, S3_msg]
# ---------------------------------------------------------------------------------------------------------



def addMetaData( metadatadir, metadata ):
    """       
      Store image03 type information into a sqlite3 database
    """
    sqlite_file = ''.join([metadatadir, '/', 'metadata.sqlite'])    # name of the sqlite database file

    # if db does not exist already we will create it
    if not os.path.isfile(sqlite_file):
            createMetaDataDB(metadatadir, metadata)
            if not os.path.isfile(sqlite_file):
                    print("Error: Could not create a database file at %s" % sqlite_file)
                    return

    table_name = 'image03'

    # Connecting to the database file
    conn = 0
    try:
        conn = sqlite3.connect(sqlite_file)
    except sqlite3.Error:
        print("Warning: Could not connect to database file %s... wait and try again." % sqlite_file)
        time.sleep(1)
        try:
                conn = sqlite3.connect(sqlite_file)
        except sqlite3.Error:
                print("Error: Could not connect to database file %s" % sqlite_file)
                return
        pass

    c = conn.cursor()

    # A) Inserts an ID with a specific value in a second column
    keys = []
    values = []
    for key in metadata:
        keys.append(key)
        if not isinstance(metadata[key], str):
            values.append(''.join(['"', str(metadata[key]), '"']))
        else:
            values.append(''.join(['"', metadata[key], '"']))


    c.execute("INSERT INTO {tn} ({cn}) VALUES ({val})".format(tn=table_name, cn=(','.join(keys)), val=(','.join(values))))

    conn.commit()
    conn.close()

def exportMetaData( filename ):
    """
       write the data from the sql table to disk as a csv file
    """
    sqlite_file = ''.join([metadatadir, '/', 'metadata.sqlite'])    # name of the sqlite database file
    
    # Connecting to the database file
    conn = sqlite3.connect(sqlite_file)
    c = conn.cursor()

    c.execute('SELECT * FROM {tn}'.format(tn="image03"))
    header = list(map(lambda x: x[0], c.description))
    #del header[0]
    all_rows = c.fetchall()
    #print(all_rows)
    with open(filename, "w") as f:
        writer = csv.writer(f)
        writer.writerow(['image','3'])
        writer.writerow(header)
        writer.writerows(all_rows)
    print('Export to %s done.' % filename)

    conn.commit()
    conn.close()

# return a value from REDCap
def getValueFromREDCap( pGUID, event, measure ):
    if not os.path.exists('config.json'):
            print("Error: Could not find local config.json, no information on how to connect to REDCap")
            sys.exit()
    token = ""
    redcap_url = ""
    with open('config.json','r') as f:
            try:
                    config = json.load(f)
            except ValueError:
                    print("Error: could not read config.json")
                    sys.exit()
            try:
                    token = config["REDCAP"][0]['token']
                    redcap_url = config["REDCAP"][0]['server_url']
            except KeyError:
                    print("Error: expected keys ['REDCAP'][0]['token'] and ['REDCAP'][0]['server_url']")
                    sys.exit()
    
    buf = io.BytesIO() # StringIO()
    # variables we need from REDCap
    data = {
        'token': token,
        'content': 'record',
        'format': 'json',
        'type': 'flat',
        'records[0]': pGUID,
        'fields[0]': measure,
        'fields[1]': 'id_redcap',
        'fields[1]': 'enroll_total',
        'events[0]': event,
        'rawOrLabel': 'raw',
        'rawOrLabelHeaders': 'raw',
        'exportCheckboxLabel': 'false',
        'exportSurveyFields': 'false',
        'exportDataAccessGroups': 'true',
        'returnFormat': 'json'
    }
    ch = pycurl.Curl()
    ch.setopt(ch.URL, redcap_url)
    ch.setopt(ch.HTTPPOST, list(data.items()))
    ch.setopt(ch.WRITEFUNCTION, buf.write)
    ch.perform()
    ch.close()
    vv = StringIO(buf.getvalue().decode('UTF-8'))
    v = json.load( vv ) # buf.getvalue() )
    buf.close()
    if isinstance(v,list) and (len(v) > 0) and ('enroll_total___1' in v[0]):
        if v[0]['enroll_total___1'] != "1":
                print("Warning: participant %s is not enrolled" % pGUID)
                return -1
    else:
        print("Warning: Could not read enroll_total from REDCap, don't know if this participant is enrolled")
        return -1
    if isinstance(v,list) and (len(v) > 0) and (measure in v[0]):
        return v[0][measure]
    return -1

# stub, not used for now
def checkHandedness( demo ):
    # check first if handedness is already part of the demo data 
    if (len(demo) > 0) and ('handedness' in demo[0]):
       return demo #nothing else to do
    print("We need to get handedness from REDCap for everyone...")
    buf = io.BytesIO() # StringIO()
    # variables we need from REDCap
    data = {
        'token': '',
        'content': 'record',
        'format': 'json',
        'type': 'flat',
        'fields[0]': 'ehi_y_ss_score',
        'fields[1]': 'id_redcap',
        'fields[2]': 'redcap_event_name',
        'events[0]': 'baseline_year_1_arm_1',
        'rawOrLabel': 'raw',
        'rawOrLabelHeaders': 'raw',
        'exportCheckboxLabel': 'false',
        'exportSurveyFields': 'false',
        'exportDataAccessGroups': 'true',
        'returnFormat': 'json'
    }
    ch = pycurl.Curl()
    ch.setopt(ch.URL, 'REDCAP URL')
    ch.setopt(ch.HTTPPOST, list(data.items()))
    ch.setopt(ch.WRITEFUNCTION, buf.write)
    ch.perform()
    ch.close()
    vv = StringIO(buf.getvalue().decode('UTF-8'))
    v = json.load( vv ) # buf.getvalue() )
    buf.close()

    # merge the data into demo
    print( "query for %s" % v )
    for entry in demo:
            for entry2 in v:
                    if entry['pGUID'] == entry2['id_redcap']:
                            if entry2['ehi_y_ss_score'] == '':
                                print("Error: could not get handedness for pGUID %s" % demo[entry]['pGUID'])
                                continue
                            if entry2['ehi_y_ss_score'] == '0':
                                entry['handedness'] = 'L'
                            elif entry2['ehi_y_ss_score'] == '1':
                                entry['handedness'] = 'R'
                            else:
                                entry['handedness'] = 'M'
    # save the demo again as a new file
    print(json.dumps(demo))

    return demo


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    ss = '\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix)
    print(ss.encode("utf-8"), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

def lookupEdat2(si_uid, eprime):
        path = ''
        eprime_ext = ''
        for i in range(0,eprime.shape[0]):
          if eprime['SeriesInstanceUID'][i] == si_uid:
            path = eprime['eprime_file_name'][i]	
            break
        if pd.isnull(path):
           return -1
        else: 
          if path.endswith('.csv'):
            print('csv e-prime file found')
            eprime_ext = 'csv'
          elif path.endswith('.txt'):
            print('txt e-prime file found')
            eprime_ext = 'txt'
          else:
            print('Error: file extension for E-Prime neither csv nor txt')
            return -1

        if not os.path.isfile(path):
            print("Error: E-Prime file could not be found - its in the spreadsheet but not on disk")
            return -1
        return (path, eprime_ext)

def bidsFolderByClassifyType( type ):
        ret = 'fmap'
        if 'ABCD-DTI' in type:
            ret = 'dwi'                  
        if ('ABCD-SST-fMRI' in type) or ('ABCD-MID-fMRI' in type) or ('ABCD-nBack-fMRI' in type) or ('ABCD-rsfMRI' in type): 
            ret = 'func'               
        if ('ABCD-T1' in type) or ('ABCD-T1-NORM' in type) or ('ABCD-T2' in type) or ('ABCD-T2-NORM' in type):
            ret = 'anat'
        return ret

#  Hauke,    May 2016               
if __name__ == "__main__":
        lfn = ''.join([ os.path.dirname(os.path.abspath(__file__)), os.path.sep, '/anonymizer.log' ])
        # logging.basicConfig(filename=lfn,format='%(levelname)s:%(asctime)s: %(message)s',level=logging.DEBUG)
        log = logging.getLogger('MyLogger')
        log.setLevel(logging.DEBUG)
        handler = logging.handlers.RotatingFileHandler( lfn, maxBytes=1e+7, backupCount=5 )
        handler.setFormatter(logging.Formatter('%(levelname)s:%(asctime)s: %(message)s'))
        log.addHandler(handler)

        inputfile = ''
        demographics = ''
        outputdir = '.'
        metadatadir = '.'
        exportfile = '.'
        eprimefile = ''
        try:
                opts,args = getopt.getopt(sys.argv[1:],"hi:d:o:m:p:e",["input=","demographics=","outputdir=","metadatadir=","eprime=","export="])
        except getopt.GetoptError:
                print('anonymizer.sh -i <input tgz> -d <demographics information in json format> -o <output directory> -m <meta data directory> -p <eprime spreadsheet>')
                sys.exit(2)
        for opt, arg in opts:
                if opt == '-h':
                        print('anonymizer.sh -i <input tgz> -d <demographics information in json format> -m <image03 meta data directory>')
                        sys.exit()
                elif opt in ("-i", "--input"):
                        inputfile = arg
                elif opt in ("-d", "--demographics"):
                        demographics = arg
                        #print(demographics)
                elif opt in ("-o", "--outputdir"):
                        outputdir = arg
                elif opt in ("-m", "--metadatadir"):
                        metadatadir = arg
                elif opt in ("-e", "--export"):
                        exportfile = arg
                        exportMetaData( exportfile )
                        os.sys.exit(0)
                elif opt in ("-p", "--eprime"):
                        eprimefile = arg  

        outputdir = os.path.abspath(outputdir)

        if demographics == '':
                print('Error: we require a demographics json file with gender, pGUID, dob')
                print('anonymizer.sh -i <input tgz> -d <demographics information in json format> -o <output directory> -m <meta data directory> -p <eprime spreadsheet>')
                log.error('Error: we require a demographics file with gender, pGUID, dob')
                sys.exit()

        demo = []
        if os.path.exists(demographics):
                with open(demographics,'r') as f:
                     try:
                        demo = json.load(f)
                     except ValueError:
                        print("Error: could not read demographics data from file")
                        log.error("Error: could not read demographics data from file")
                        sys.exit(0)
        if demo == []:
                print('Error: could not read demographics as json')
                log.error('Error: could not read demographics as json')
                sys.exit(0)
        
        # handedness should be in there, if not we need to compute it once
        # demo = checkHandedness( demo )

        if not os.path.exists(outputdir):
                print("Warning: output directory %s does not exist, try to create..." % outputdir)
                os.mkdir(outputdir)
                if not os.path.exists(outputdir):
                        print("Error: could not create output directory %s" % outputdir)
                        log.error("Error: could not create output directory %s" % outputdir)
                        sys.exit(0)

        if not os.path.exists(metadatadir):
                print("Warning: metadata directory %s does not exist, try to create..." % outputdir)
                os.mkdir(metadatadir)
                if not os.path.exists(metadatadir):
                        print("Error: could not create metadata directory %s" % metadatadir)
                        log.error("Error: could not create metadata directory %s" % metadatadir)
                        sys.exit(0)
        
        if eprimefile == '':
                print('Error: we require an e-prime file')
                print('anonymizer.sh -i <input tgz> -d <demographics information in json format> -o <output directory> -m <meta data directory> -p <eprime spreadsheet>')
                log.error('Error: we require an e-prime file')
                sys.exit()

        eprime = []
        if os.path.exists(eprimefile):
                with open(eprimefile,'r') as f:
                     try:
                        eprime = pd.read_csv(f)
                     except ValueError:
                        print("Error: could not read eprime data from file")
                        log.error("Error: could not read eprime data from file")
                        sys.exit(0)
        
        if eprime.shape[0] == 0:
                print('Error: could not read e-prime file as csv')
                log.error('Error: could not read e-prime file as csv')
                sys.exit(0)


        if 1:
                # general lookupEdat2purpose disclaimer
                if not os.path.exists(inputfile):
                    print("Error: input argument is neither directory nor file")
                    log.error("Error: input argument is neither directory nor file")
                    sys.exit(0)

                # get the values that we will use for anonymizing 
                anonInfo = {}
                anonInfo['inputfile'] = inputfile
                # get the subject ID and event name
                # we have a pattern for the subject ID (only accept NDAR pGUIDS)
                pattern  = re.compile('(NDAR_INV[0-9A-Z]+)_(.*(?=_Session))')
                pattern2 = re.compile('(ABCD[0-9A-Za-z]+)_(.*(?=_Session))')
                pattern3 = re.compile('(PhantomTravelingHuman_[0-9]+_[A-Z]+)_(.*(?=_Session))')
                pattern4 = re.compile('(TestQA_[A-Z]+)_(.*(?=_Session))')
                pattern5 = re.compile('.*TEST.*')
                fn, fext = os.path.splitext(os.path.basename(inputfile))
                #print("test: %s" % fn)
                gg = re.search(pattern, fn)
                if gg == None:
                    gg = re.search(pattern2, fn)
                if gg == None:
                    gg = re.search(pattern3, fn)
                if gg == None:
                    gg = re.search(pattern4, fn)
                # test for a TEST pGUID, has 'TEST' in the pGUID
                gg2 = re.search(pattern5, fn)
                if gg2 != None:
                    # if we find a TEST, don't share
                    print("Warning: appears to be a TEST pGUID, don't share")
                    log.error("Warning: appears to be a TEST pGUID %s, don't share" % (fn))
                    sys.exit(0)

                if gg == None:
                    print("Could not match the participant and event name for %s" % (fn))
                    log.error("Could not match the participant and event name for %s" % (fn))
                    sys.exit(0)
                matches = gg.groups()
                #print("found this many entries: %d from %s %s" % (len(matches), fn, matches[0]))
                if len(matches) == 2:
                    anonInfo["pGUID"] = matches[0]
                    anonInfo["event"] = matches[1]
                    anonInfo["pGUID_BIDS"] = re.sub(r'_([a-zA-Z0-9])', lambda m: m.group(1).upper(), anonInfo["pGUID"])
                    anonInfo["event_BIDS"] = re.sub(r'_([a-zA-Z0-9])', lambda m: m.group(1).upper(), anonInfo["event"])
                else:
                    print("Error: Could not match the pGUID and name of the session for %s" % (fn))
                    log.error("Error: Could not match the pGUID and name of the session for %s" %(fn))
                    sys.exit(0)

                # find out what the radiology review score is for this participant
                score = getValueFromREDCap( anonInfo['pGUID'], anonInfo['event'], 'mrif_score' )
                #print('score:', score) 
                if not ((score == "1") or (score == "2")):
                    print("Error: participant score %s for %s not in allowed range (1,2)" % (score, anonInfo['pGUID']))
                    log.error("Error: participant %s score not in allowed range (%s)" %( anonInfo['pGUID'], score))
                    sys.exit(0)

                # in order to find out what the type of scan is we need to look at the json file that comes with this tgz (should be inside)


                # we have two modes, either the input is a directory (mode = dir) or a tgz file (mode = file)
                mode = ""
                if os.path.isdir(inputfile):
                    mode = "dir"
                    print("Error: Directory mode is disabled for now, provide a TGZ for sharing.\n")
                    log.error("Error: Directory mode is disabled for now, provide a TGZ for sharing.\n")
                    sys.exit(0)
                if os.path.isfile(inputfile) and tarfile.is_tarfile(inputfile):
                    mode = "file"
                if mode == "":
                    print("Error: Input is neither tar nor directory")
                    log.error("Error: input is neither tar nor directory")
                    sys.exit(0)

                # read the tags to delete/change
                tagsFile = ''.join([ os.path.dirname(os.path.abspath(__file__)), os.path.sep, 'anonymizeTags.json'])
                anonymizeThis = []
                if os.path.exists(tagsFile):
                        with open(tagsFile,'r') as f:
                                try:
                                        anonymizeThis = json.load(f)
                                except ValueError:
                                        print("Error: could not read anonymizeTags.json file, syntax error")
                                        log.error("Error: could not read anonymizeTags.json file, syntax error")
                                        sys.exit(0)
                else:
                        log.critical("Error: could not find the list of tags to remove (anonymizeTags). This json file is expected in the local path.")
                        sys.exit(0)

                # create a pool of worker to speed things up
                # (tests show that with 6 CPU's we will still need about 2min to anonymize 14,000 files of a study,
                #  with 2 CPUs we need 4min)
                pool = Pool(processes=1)
                res  = []
                if mode == "dir":
                    for r,d,f in os.walk(inputfile):
                            for file in f:
                                    df = os.path.join(inputfile,file)
                                    #print df
                                    res.append(pool.apply_async(anonymize, args=(df,)))
                    for proc in res:
                        proc.get(timeout=1)

                if mode == "file":
                    # if we have a tgz file we should use tarfile to read it
                    # we can write the new data into a new tar?
                    fn, fext = os.path.splitext(inputfile)
                    fnjson = ''.join([ fn, '.json'])
                    # try to read the json file that describes the scan
                    tarin = tarfile.open(inputfile, 'r:*')
                    # we are getting the members of the tar file here. this is expensive - reads the whole file
                    members = tarin.getmembers()
                    reT = re.compile(r'.*[.]+hdr')
                    jfile = [m for m in members if reT.search(m.name)]
                    if len(jfile) == 1:
                            # found a hdr file, skip this package
                            print("Error: detected hdr file inside tgz, assume %s contains k-space, don't do anything" % anonInfo['inputfile'])
                            log.error("Error: detected hdr file inside %s, assume this is k-space, don't do anything" % anonInfo['inputfile'])
                            os.sys.exit(0)

                    # lets try to find the json file inside that will tell us what sort of data this is
                    reT = re.compile(r'.*[.]+json')
                    jfile = [m for m in members if reT.search(m.name)]
                    if len(jfile) == 1:
                        fobj = tarin.extractfile(jfile[0])
                        # read the json data and convert to string for json decode
                        str_fobj = fobj.read().decode("utf-8")
                        data = json.loads(str_fobj)
                        # we can check if this dataset has been anonymized before, should contain a Anonymized entry
                        if 'Anonymized' in data:
                                print("Detected already anonymized dataset (Anonymized tag present in included json). Ignoring teprime_flags = lookupEdat2(anonInfo['SeriesInstanceUID']) his tgz.")
                                log.error("Detected already anonymized dataset (Anonymized tag present in included json). Ignoring this tgz.")
                                os.sys.exit(-1)

                        anonInfo['StudyDate'] = ""
                        anonInfo['StudyTime'] = ""
                        anonInfo['SeriesDate'] = ""
                        anonInfo['SeriesTime'] = ""
                        anonInfo['Manufacturer'] = ""
                        anonInfo['ManufacturerModelName'] = ""
                        anonInfo['SeriesDescription'] = ""
                        anonInfo['SeriesInstanceUID'] = ""
                        try:
                                anonInfo['StudyDate'] = data['StudyDate']
                        except KeyError:
                                pass
                        try:
                                #eprime_flags = lookupEdat2(anonInfo['SeriesInstanceUID'], eprime)
                                #if eprime_flags == -1:
                                #        print("Error: Event e-prime file not found")
                                #        log.error("Error: Event e-prime file not found for participant %s " %( anonInfo['pGUID'] ))
                                #        sys.exit(0) 
                                anonInfo['SeriesInstanceUID'] = data['SeriesInstanceUID']
                        except KeyError:
                                pass
                        try:
                            anonInfo['SeriesInstanceUID'] = data['SeriesInstanceUID']
                        except KeyError:
                                pass
                        try:
                            anonInfo['StudyTime'] = data['StudyTime']
                        except KeyError:
                                pass
                        anonInfo['SeriesDate'] = anonInfo['StudyDate']
                        anonInfo['SeriesTime'] = anonInfo['StudyTime']
                        try:
                            anonInfo['SeriesDate'] = data['SeriesDate']
                        except KeyError:
                                pass
                        try:
                            anonInfo['SeriesTime'] = data['SeriesTime']
                        except KeyError:
                                pass
                        try:
                            anonInfo['Manufacturer'] = data['Manufacturer']
                        except KeyError:
                                pass
                        try:
                                #eprime_flags = lookupEdat2(anonInfo['SeriesInstanceUID'], eprime)
                                #if eprime_flags == -1:
                                #        print("Error: Event e-prime file not found")
                                #        log.error("Error: Event e-prime file not found for participant %s " %( anonInfo['pGUID'] ))
                                #        sys.exit(0) 

                                anonInfo['ManufacturerModelName'] = data['ManufacturerModelName']
                        except KeyError:
                                pass
                        try:
                            anonInfo['SeriesDescription'] = data['SeriesDescription']
                        except KeyError:
                                pass

                        # we should have a section in here about the ClassifyType
                        anonInfo["ABCDType"] = "unknown"
                        anonInfo["ClassifyType"] = []
                        if 'ClassifyType' in data:
                                anonInfo["ClassifyType"] = data["ClassifyType"]
                                anonInfo["ABCDType"] = "_".join(s for s in data["ClassifyType"] if "ABCD".lower() in s.lower())
                                
                                if anonInfo["ABCDType"] == "unknown":
                                        # could be QA?
                                        anonInfo["ABCDType"] = "_".join(s for s in data["ClassifyType"] if "QA".lower() in s.lower())
                        if anonInfo["ABCDType"] == "":
                                anonInfo["ABCDType"] = "unknown"

                        # now we can add the demographic information for this participant
                        anonInfo['gender'] = ''
                        anonInfo['dob'] = ''
                        for val in demo:
                                if (val['pGUID'] == anonInfo['pGUID']) and (val['dob'] != ''):
                                        anonInfo['gender'] = val['gender']
                                        anonInfo['dob'] = val['dob']
                                        break
                    else:
                            print("Error: could not find a single json file in the tgz, quit here (maybe we should keep going?)")
                            os.sys.exit(0)

                    #print('ABCDType', anonInfo['ABCDType'])
                    #print('other' , anonInfo['ManufacturerModelName'] ) 
                    
                    if (anonInfo['ABCDType'] == 'unknown') or (anonInfo['ABCDType'] == 'ABCD-Physio'):
                        # bail out, don't share
                        print("Error: ClassifyType not shared %s" % (" + ".join(anonInfo["ClassifyType"])))
                        log.error("Error: ClassifyType not shared %s" % (" + ".join(anonInfo["ClassifyType"])))
                        os.sys.exit(0)
                    
                    #calculate a series time without the fractional seconds
                    seriestime = anonInfo['SeriesTime']
                    if "." in seriestime:
                            seriestime = seriestime.split(".")[0] 

                    # find out if the eprime event file exist 
                    # check this only for the task runs 
                    eventEdatFileName_description = '' 
                    if 'ABCD-SST-fMRI' in anonInfo['ClassifyType'] or 'ABCD-MID-fMRI' in anonInfo['ClassifyType'] or 'ABCD-nBack-fMRI' in anonInfo['ClassifyType']: 
                        eprime_flags = lookupEdat2(anonInfo['SeriesInstanceUID'], eprime) 
                        if eprime_flags == -1:
                                print("Error: Event e-prime file not found")
                                log.error("Error: Event e-prime file not found for participant %s " %( anonInfo['pGUID'] ))
                                sys.exit(0)    
                        else: 
                                bidstype = bidsFolderByClassifyType( anonInfo['ClassifyType'] )
                                eventEdatFileName = "sub-%s/ses-%s/%s/%s_run-%s%s-EventRelatedInformation.%s" % (anonInfo['pGUID_BIDS'],
                                                                                                                 anonInfo['event_BIDS'],
                                                                                                                 bidstype,
                                                                                                                 anonInfo['ABCDType'],
                                                                                                                 anonInfo['SeriesDate'],
                                                                                                                 seriestime,
                                                                                                                 eprime_flags[1])
                                eventEdatFileName_description = 'subject-level task response'                                                      

                    #study date NIH format              
                    sday2 = datetime.datetime.strptime(anonInfo['StudyDate'],'%Y%m%d').strftime('%m/%d/%Y')

                    # our record to update contains {"name": "interview_date", "value": "04/06/2017"}
                    # but NDA requires {"name": "interview_date", "value": "04/06/2017 00:00:00"}, so:
                    if len(anonInfo['StudyTime']) >= 6:
                        sday2 = sday2 + ' ' + anonInfo['StudyTime'][0:2] + ':' + anonInfo['StudyTime'][2:4] + ':' + anonInfo['StudyTime'][4:6]
                    else:
                        sday2 = sday2 + " 00:00:00"

                    outtarname = ''.join([ os.path.abspath(outputdir), os.path.sep, anonInfo['pGUID_BIDS'], '_', 
                                           anonInfo['event_BIDS'], '_', anonInfo["ABCDType"], '_', anonInfo['SeriesDate'], 
                                           seriestime, '.tgz'])
                    print("Write to %s ..." % outtarname)
                    log.info("Write to %s ..." % outtarname)

                    dti_flag = ''
                    scan_type = 'Field Map'
                    if 'ABCD-DTI' in anonInfo['ClassifyType']:
                        dti_flag = 'Yes'
                        scan_type = 'multi-shell DTI'
                    if ('ABCD-T1' in anonInfo['ClassifyType']) or ('ABCD-T1-NORM' in anonInfo['ClassifyType']):
                        scan_type = 'MR structural (T1)'        
                    if ('ABCD-T2' in anonInfo['ClassifyType']) or ('ABCD-T2-NORM' in anonInfo['ClassifyType']):
                        scan_type = 'MR structural (T2)'        


                    id_flag = ''
                    if 'ABCD-SST-fMRI' in anonInfo['ClassifyType']: 
                        id_flag = '650'
                        scan_type = 'fMRI'
                    if 'ABCD-MID-fMRI' in anonInfo['ClassifyType']:
                        id_flag = '648'
                        scan_type = 'fMRI'
                    if 'ABCD-nBack-fMRI' in anonInfo['ClassifyType']: 
                        id_flag = '651'
                        scan_type = 'fMRI'
                    if 'ABCD-rsfMRI' in anonInfo['ClassifyType']: 
                        id_flag = '649'                 
                        scan_type = 'fMRI'

                    if os.path.exists(outtarname):
                            print("Error: file already exists, please move away before re-compressing...")
                            log.error("output file already exists, stop processing.")
                            os.sys.exit(0)
                    tarout = tarfile.open(outtarname, 'w:gz')
                    i = 0
                    for entry in tarin.getmembers():
                        en = tarin.extractfile(entry)
                        if hasattr(en, 'read'):
                          fobj = en.read()
                          if fobj == None:
                                  continue
                          anonInfo['origFileName'] = entry
                          # save the file temporarily (which is stupid, should instead work with a in memory file such as StringIO)
                          with tempfile.NamedTemporaryFile() as temp:
                                  try:
                                          temp.write(fobj)
                                  except OSError:
                                          pass
                                  temp.flush()
                                  printProgressBar(i,len(members), prefix='Progress:', suffix=('(%d/%d) Complete' % (i,len(members)-1)), length = 50)
                                  i = i + 1
                                  # Some fields will be added to anonInfo by this call
                                  anonymize(temp.name, anonInfo, "file", tarout)
                          #memory_file = DicomBytesIO(fobj)
                          #print("is this file little endian?", memory_file)
                          #anonymize(memory_file,"file",tarout)
                          #res.append(pool.apply_async(anonymize, args=(memory_file,"file",tarout,)))
                    #for proc in res:
                    #    proc.get(timeout=1)
                    tarin.close()
                    tarout.close()

                    # Assembly meta-data record to be savev to NDA and our local database
                    new_record = {'subjectkey': anonInfo['pGUID'], # required
                                  'src_subject_id': anonInfo['pGUID'], # required
                                  'interview_date': sday2, # required
                                  'interview_age': round(float(anonInfo['PatientsAge'])*12), # required
                                  'gender': anonInfo['gender'], # required
                                  'comments_misc': anonInfo['SeriesDescription'],
                                  'image_file': outtarname, # required
                                  'image_thumbnail_file': '',
                                  'image_description': anonInfo['ABCDType'], # required DTI, fMRI, Fast SPGR, phantom
                                  'experiment_id': id_flag, # required if fMRI
                                  'scan_type': scan_type,  # required: MR diffusion; fMRI; MR structural (MPRAGE); MR structural (T1); MR structural (PD); MR structural (FSPGR); MR structural (T2); PET; ASL; microscopy; MR structural (PD, T2); MR structural (B0 map); MR structural (B1 map); single-shell DTI; multi-shell DTI; Field Map; X-Ray
                                  'scan_object': "Live", # required "Live", "Phantom"
                                  'image_file_format': "DICOM", # required
                                  'data_file2': '', #eventEdatFileName,
                                  'data_file2_type': '', #eventEdatFileName_description,
                                  'image_modality': "MRI", #required
                                  'scanner_manufacturer_pd': anonInfo['Manufacturer'], # required  
                                  'scanner_type_pd': anonInfo['ManufacturerModelName'], # required 
                                  'scanner_software_versions_pd': anonInfo['SoftwareVersion'], # required  
                                  'magnetic_field_strength': '3', 
                                  'mri_repetition_time_pd': '%.3f' % (float(anonInfo['RepetitionTime'])/1000.0), # required
                                  'mri_echo_time_pd': '%.3f' % (float(anonInfo['EchoTime'])/1000.0), # required 
                                  'flip_angle': anonInfo['FlipAngle'], # required 
                                  'acquisition_matrix': anonInfo['AcquisitionMatrix'], # required  
                                  'mri_field_of_view_pd': anonInfo['FoV'], # required 
                                  'patient_position': anonInfo['PatientPosition'], # required 
                                  'photomet_interpret': anonInfo['PhotometricInterpretation'], # required
                                  'receive_coil': anonInfo['ReceiveCoilName'],
                                  'transmit_coil': anonInfo['TransmitCoilName'],  
                                  'transformation_performed': "No", #required
                                  'transformation_type': '',
                                  'image_history': '',
                                  'image_num_dimensions': '',
                                  'image_extent1': '', 
                                  'image_extent2': '',
                                  'image_extent3': '',
                                  'image_extent4': '',
                                  'extent4_type': '',
                                  'image_extent5': '',
                                  'extent5_type': '',
                                  'image_unit1': '',
                                  'image_unit2': '',
                                  'image_unit3': '',
                                  'image_unit4': '',
                                  'image_unit5': '',
                                  'image_resolution1': '', 
                                  'image_resolution2': '',
                                  'image_resolution3': '',
                                  'image_resolution4': '',
                                  'image_resolution5': '',
                                  'image_slice_thickness': anonInfo['SliceThickness'],
                                  'image_orientation': '',
                                  'qc_outcome': '',
                                  'qc_description': '',
                                  'qc_fail_quest_reason': '',
                                  'decay_correction': '',
                                  'frame_end_times': '',
                                  'frame_end_unit': '',
                                  'frame_start_times': '',
                                  'frame_start_unit': '',
                                  'pet_isotope': '',
                                  'pet_tracer': '',
                                  'time_diff_inject_to_image': '',
                                  'time_diff_units': '',
                                  'pulse_seq': '',
                                  'slice_acquisition': '',
                                  'software_preproc': '',
                                  'study': '',
                                  'week': '',
                                  'experiment_description': '',
                                  'visit': anonInfo['event'],
                                  'slice_timing': '',
                                  'bvek_bval_files': dti_flag,
                                  'bvecfile': '',
                                  'bvalfile': '',
                                  'deviceserialnumber': '',
                                  'procdate': ''}

                # Write an entry to our meta-data file
                [miNDA_ok, miNDA_msg, S3_ok, S3_msg] = uploadToNDA( metadatadir, new_record )
                
                print('\n[miNDA_ok =', miNDA_ok)
                print('miNDA_msg: ', miNDA_msg, '\n')
                print('S3_ok =', S3_ok)
                print('S3_msg: ', S3_msg)
                if len(S3_msg) == 0:
                    S3_msg = ''
                print('S3_msg: ', S3_msg)
                print('\n' )

                new_record['miNDA_ok']  = miNDA_ok
                new_record['miNDA_msg'] = miNDA_msg
                new_record['miNDA_msg'] = new_record['miNDA_msg'].replace('\"','')
                new_record['miNDA_msg'] = new_record['miNDA_msg'].replace('\'','')

                new_record['S3_ok']  = S3_ok
                new_record['S3_msg'] = S3_msg
                new_record['S3_msg'] = new_record['S3_msg'].replace('\"','')
                new_record['S3_msg'] = new_record['S3_msg'].replace('\'','')

                addMetaData( metadatadir, new_record )

                pool.close()
                sys.exit(0)
        else:
                print("Anonymize a tgz files with DICOMs and a single .json file.")
                print("Requires as an argument the tgz file to analyze.")
                sys.exit(2)
