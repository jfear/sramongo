#!/bin/bash
################################################################################
# Uses the NCBI Entrez Direct Utilities to download results from a query to the
# SRA database. If using the conda environment the entres-direct should already


if [ ! -z "${1+x}" ]; then
    # Some of Lee's data
    efetch -db sra -id SRR3001915 -format xml | xmllint -format - > $1/sra_SRR3001915.xml
    efetch -db sra -id SRX1483046 -format xml | xmllint -format - > $1/sra_SRX1483046.xml
    efetch -db biosample -id SAMN02981965 -format xml | xmllint -format - > $1/biosample_SAMN02981965.xml
    efetch -db bioproject -id PRJNA258012 -format xml | xmllint -format - > $1/bioproject_PRJNA258012.xml

    # Some of my data
    efetch -db sra -id SRR1945105 -format xml | xmllint -format - > $1/sra_SRR1945105.xml

    # Non GEO data
    efetch -db sra -id SRR5003988 -format xml | xmllint -format - > $1/sra_SRR5003988.xml
else
    echo "Please provide the directory where to save test data."
fi
