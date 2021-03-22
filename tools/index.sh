#!/usr/bin/env bash
$ref=$1
$tool_name=$2
soft=./softwares.config
tool=`grep $tool_name"=" $soft |cut -d '=' -f 2`
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\n":$0.upper() }' $ref > ref.fa
####statistics####
genome_size=`awk '$0 !~ /^>/{print length}' ref.fa |awk '{s+=$1}END{print s}'`
A=`grep -v '>' ref.fa |grep -o 'A' | wc -l`
T=`grep -v '>' ref.fa |grep -o 'T' | wc -l`
C=`grep -v '>' ref.fa |grep -o 'C' | wc -l`
G=`grep -v '>' ref.fa |grep -o 'G' | wc -l`
N=`grep -v '>' ref.fa |grep -o 'N' | wc -l`
GC=${G}+${C}
atcg=${A}+${T}+${C}+${G}
GC_content=`echo "scale=2; ${GC}*100/${atcg}" |bc`
echo -e "A:${A}\nT:${T}\nC:${C}\nG:${G}\nGC content(No N's):${GC_content}%\nTotal:{genome_size}" > ref.fa.stat.txt
####index####
if [[ $tool_name = 'BWA' && $genome_size -lt 4000000000 ]] ; then  ## <4G small genome
  $tool index ref.fa -p ref
elif [[ $tool_name = 'BWA' && $genome_size -ge 4000000000 ]] ; then ## >4G large genome
  $tool index -a bwtsw -p ref
elif [[ $tool_name = 'samtools' ]] ;then
  $tool faidx ref.fa
elif [[ $tool_name = 'hisat2' ]];then
  $tool"-build" -p 8 ref.fa ref
elif [[ $tool_name = 'picard' ]];then
  $tool  CreateSequenceDictionary R=ref.fa O=ref.dict
else
  echo "Error,no such tool in softwares.config or input wrong name, please add software in software.config or use correct name"
fi


    



