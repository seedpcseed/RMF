#!/bin/sh

if [ -d "log" ]; then
   mkdir log
fi

echo "#---- Starting up RMetaflow"
echo ""
echo "#---- Have a great day!"
echo ""
sbatch RMetaflow.sh 