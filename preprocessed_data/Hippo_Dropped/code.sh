#!/bin/bash

bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/rnaseq-run-all.sh --experiment "HIPPO_riboZero_BrainSeq_Phase2" --prefix "droppedSamples" --reference "hg38" --stranded "reverse" --ercc "TRUE" 
bash /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step8-callVariants.sh --experiment "HIPPO_riboZero_BrainSeq_Phase2" --prefix "droppedSamples" --reference "hg38" 