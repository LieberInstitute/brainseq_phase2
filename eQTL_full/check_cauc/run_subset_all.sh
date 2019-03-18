#!/bin/bash

mkdir -p logs

qsub run_subset_dlpfc.sh
qsub run_subset_hippo.sh
qsub run_subset_interaction.sh
