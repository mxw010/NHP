#!/bin/bash

for f in `ls ../raw_data/*.fastq.gz`; do
    f=`basename $f`
    output=`echo $f | sed 's/.fastq.gz//g'`
    ../test/fastq_nnk_extractor \
    "gtcctatggacaagtggccacaaaccaccagagtgcccaaNNKNNKNNKNNKNNKNNKNNKgcacaggcgcagaccggctgggttcaaaaccaaggaatacttcc" \
    ../raw_data/${f} \
    0 \
    ${output}.txt \
    10 >> summary.txt
done


cut -d":" -f2 QC.txt | tr "\n" " " | sed 's/_001.fastq.gz//g' | sed 's/..\/raw_data//g' | tr "/" "\n" | sed '/^$/d' | tr "_" " " | tr "-" " " | cut -d" " -f2,7- > temp
echo "File Total_Reads QC Matching_Flanking_Sequences Valid_NNK Invalid_NNK Final_Extracted" > QC.txt
cat temp >> QC.txt
sed '/^$/d' QC.txt