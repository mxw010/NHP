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
for (( i=1; i<=`wc -l ${work_dir}/data/amplicons.txt | awk '{print $1}'`; i++ )); do
    amplicon=`head -$i ${work_dir}/data//amplicons.txt | tail -1 | awk '{print $1}'`
    seq=`head -$i ${work_dir}/data//amplicons.txt | tail -1 | awk '{print $2}'`
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

# basic QC processing

mkdir processed
for sample in `ls *_output.txt | sed 's/_R[12]_output.txt//g' | sort | uniq`; do
    Rscript /igm/home/mxw010/NHP/scripts/QC_filtering.r $sample processed/${sample}_clean.txt
done

## AAV6: forward 75, reverse 33
## AAV9: forward 40, reverse 115
## RH10: forward 85, reverse 42

cd processed 

# numerical summary of post QC reads (invalid NNK and stop codons)
for file in `ls *_clean.txt`; do
    capsid=`echo $file | cut -f1 -d"_"`
    sample=`echo $file | cut -f2 -d"-"`
    echo -n "$capsid $sample " >> summary.csv
    awk  -F',' '{n++
        if ($7 == "FALSE") invalid++
    if ($8 == "TRUE") stop++
    }
    END {
    print n, invalid,  stop
    }' $file >> summary.csv
done

tr " " "," < summary.csv > temp
mv temp summary.csv

# get unique peptides count
echo "Sample,Capsid,Pepetide_Counts" > peptide_counts.csv
for file in `ls *.txt`; do
    count=`cut -f4 -d"," $file | grep  -v "\*" | sort | uniq | wc -l`
    count=`echo $count-1 | bc -l`
    capsid=`echo $file | cut -f1 -d"_"`
    sample=`echo $file | cut -f2 -d"-"`
    echo "$sample,$capsid,$count" >> peptide_counts.csv
done