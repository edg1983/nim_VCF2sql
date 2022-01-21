#!/bin/bash
# author Edoardo Giacopuzzi
# template single job submission

#$ -q short.qc
#$ -cwd -j y
#$ -o cohort_deepvariant_bulk_null.log
#$ -N sqllite
#$ -l h_vmem=4G 
#$ -pe ramdisk 5

echo "started on `hostname` at `date`"

./sqllite_tools_bulk \
--input cohort_deepvariant2.v2r.idx.tsv.gz \
--outdb ${TMPDIR}/cohort_deepvariant_bulk_null.db \
-v -g

mv ${TMPDIR}/cohort_deepvariant_bulk_null.db ./

echo "Stop " $(date)
