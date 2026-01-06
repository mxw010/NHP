#!/bin/bash


# compile C++
g++ -O3 -o build/fastq_nnk_extractor src/NNK_sequence_extractor.cpp -lz -std=c++11

# for debugging
# g++ -Wall -o fastq_nnk_extractor fastq_nnk_extractor.cpp -lz -std=c++11



g++ -Wall -o fastq_nnk_extractor ../src/NNK_sequence_extractor.cpp -lz -std=c++11

./fastq_nnk_extractor "gtcctatggacaagtggccacaaaccaccagagtgcccaaNNKNNKNNKNNKNNKNNKNNKgcacaggc" temp.txt.gz 0 output.txt 10

## Sciatic Nerve
work_dir=`pwd`
cd data
mkdir -p sciatic_nerve
cd sciatic_nerve
for (( i=1; i<=`wc -l ../amplicons.txt`; i++)); do
    amplicon=`head -$i ../amplicons.txt | tail -1 | awk '{print $1}'`
    seq=`head -$i ../amplicons.txt | tail -1 | awk '{print $2}'`
    for f in `ls ${work_dir}/fastq/sciatic_nerve/`; do
        f=`basename $f`
        output=`echo $f | sed 's/_001.fastq.gz//g'`
        ${work_dir}/build/fastq_nnk_extractor \
        $seq \
        ${work_dir}/fastq/sciatic_nerve/${f} \
        0 \
        ${work_dir}/data/sciatic_nerve/${amplicon}_${output}_output.txt \
        10 >> ${work_dir}/data/sciatic_nerve/${amplicon}_${output}_summary.txt
    done
done


cut -d":" -f2 QC.txt | tr "\n" " " | sed 's/_001.fastq.gz//g' | sed 's/..\/raw_data//g' | tr "/" "\n" | sed '/^$/d' | tr "_" " " | tr "-" " " | cut -d" " -f2,7- > temp
echo "File Total_Reads QC Matching_Flanking_Sequences Valid_NNK Invalid_NNK Final_Extracted" > QC.txt
cat temp >> QC.txt
sed '/^$/d' QC.txt