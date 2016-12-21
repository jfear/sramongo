#!/bin/bash
################################################################################
# Uses the NCBI Entrez Direct Utilities to download results from a query to the
# SRA database. If using the conda environment the entres-direct should already


if [ ! -z "${1+x}" ]; then
    # Some of Lee's data
    efetch -db sra -id SRR3001915 -format xml | xmllint -format - > $1/sra_SRR3001915.xml
    efetch -db biosample -id SAMN02981965 -format xml | xmllint -format - > $1/biosample_SAMN02981965.xml
    efetch -db bioproject -id PRJNA258012 -format xml | xmllint -format - > $1/bioproject_PRJNA258012.xml
    efetch -db pubmed -id 26732976 -format xml | xmllint -format - > $1/pubmed_26732976.xml

    # Some of my data
    efetch -db sra -id SRR1945105 -format xml | xmllint -format - > $1/sra_SRR1945105.xml

    # New data
    efetch -db sra -id SRR5100239 -format xml | xmllint -format - > $1/sra_SRR5100239.xml
    efetch -db sra -id ERR1662611 -format xml | xmllint -format - > $1/sra_ERR1662611.xml

    # modENCODE Example
    efetch -db sra -id SRR1197474 -format xml | xmllint -format - > $1/sra_SRR1197474_mod.xml
else
    echo "Please provide the directory where to save test data."
fi
