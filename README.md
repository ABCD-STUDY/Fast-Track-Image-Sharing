# DICOM Data Image Sharing

The Adolescent Brain Cognitive Development study (ABCD) is sharing raw DICOM images in a BIDS like file directory format. This script will anonymize DICOM files, add meta-data and create a TGZ file suitable for sharing on a platform such as NDA (National Data Archive).

The Brain Imaging Data Structure (http://bids.neuroimaging.io) describes a file format for sharing of nifti-format encoded volumetric data. ABCD is using its directory structure to share DICOM data. In a trivial process (dcm2nii) the data generated by this tool can be converted to a fully BIDS compatible structure. As is BIDS tools might not be able to process the provided data.

## Input Data

Data exists on the local file system packaged as TGZ encoded directory that contain all DICOM files of a single image series. Together with the DICOM files the TGZ also contains a single json file (extension .json) that describes the data. This file is created by the processSingleFile python script available on the FIONASITE repository (FIONASITE/server/bin/processSingleFile.py). Here an example for the content of the generated file:
```json
{
    "AccessionNumber": "",
    "AcquisitionLength": "TA 06:00",
    "AcquisitionMatrix": "[90, 0, 0, 90]",
    "ActiveCoils": "HEA;HEP",
    "Anonymized": "fast track v1.0, @DAIC April 2017",
    "ClassifyType": [
        "SIEMENS",
        "mosaic",
        "original",
        "ABCD-SST-fMRI"
    ],
    "DateOfBirth": "20020215",
    "EchoTime": "30",
    "ImageType": "['ORIGINAL', 'PRIMARY', 'M', 'ND', 'MOSAIC']",
    "IncomingConnection": "removed",
    "InstanceNumber": "445",
    "Manufacturer": "SIEMENS",
    "ManufacturerModelName": "Prisma_fit",
    "Modality": "MR",
    "NumFiles": "removed",
    "PatientID": "NDAR_INVXXXXXXXX",
    "PatientsAge": "116",
    "PatientsSex": "M",
    "PhaseEncodingDirectionPositive": "0",
    "RepetitionTime": "800",
    "SOPClassUID": "MR Image Storage",
    "ScanningSequence": "EP",
    "SequenceName": "epfSM2d1_90",
    "SequenceVariant": "SK",
    "SeriesDescription": "SIEMENS, mosaic, original, ABCD-SST-fMRI",
    "SeriesInstanceUID": "XXX",
    "SeriesNumber": "28",
    "SeriesTime": "122945.581000",
    "SliceLocation": "-59.027603148797",
    "SliceSpacing": "2.3999999741376",
    "SliceThickness": "2.4000000953674",
    "StudyDate": "20170122",
    "StudyDescription": "Adolescent Brain Cognitive Development Study",
    "StudyInstanceUID": "XXX",
    "StudyTime": "112752.432000",
    "siemensDiffusionInformation": [],
    "siemensUUID": "X-X-X"
}
```

## Output

The script will create an output TGZ file that contains a directory structure that follows BIDS closely.

```
sub-NDARINVXXXXXXXX/
  ses-baselineYear1Arm1/
    func/
      ABCD-SST-FMRI_run-20000000/
      ABCD-SST-fMRI_run-20000000-EventRelatedInformation.txt
      ABCD-SST-fMRI_run-20000000.json
```

## Workflow

Together with the output TGZ the script will store NDA's image03 information in a sqlite database. After anonymizing a sufficient number of image series a call of anonymizer.sh with the option '-e' will export this information as a csv file suitable for the NDA upload tool.
