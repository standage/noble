#!/usr/bin/env bash
set -eo pipefail
shopt -s expand_aliases

alias vcf2bed="awk -v OFS='\t' '{ print \$1, \$2-2, \$2+2 }'"

refr=$1
pred=$2

correct=$(
    bedtools intersect \
            -a <(grep ',0/0,0/0' $refr | vcf2bed) \
            -b <(vcf2bed $pred) \
            -u -wa \
        | wc -l
)
missing=$(
    bedtools intersect \
            -a <(grep ',0/0,0/0' $refr | vcf2bed) \
            -b <(vcf2bed $pred) \
            -v -wa \
        | wc -l
)
wrong=$(
    bedtools intersect \
            -b <(grep ',0/0,0/0' $refr | vcf2bed) \
            -a <(vcf2bed $pred) \
            -v -wa \
        | wc -l
)

echo "Correct: $correct"
echo "Missing: $missing"
echo "Wrong:   $wrong"
