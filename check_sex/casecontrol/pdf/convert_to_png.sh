#!/bin/bash

#https://askubuntu.com/questions/50170/how-to-convert-pdf-to-image
## Ran this on my laptop
for i in *pdf;
    do
    echo ${i}
    # https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
    NAME=`echo "${i}" | cut -d'.' -f1`
    echo ${NAME}
    convert -density 300 ${i} -quality 100 ${NAME}.png
done

## Upload back to JHPCE
scp *png e:/dcl01/lieber/ajaffe/lab/brainseq_phase2/check_sex/casecontrol/pdf/
