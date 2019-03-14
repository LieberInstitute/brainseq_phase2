#!/bin/bash

# Usage:
# bash run_all.sh

mkdir -p eqtl_tables
mkdir -p eqtl_tables/logs
mkdir -p rdas

qsub eQTL_dlpfc.sh
qsub eQTL_hippo.sh
qsub eQTL_interaction.sh