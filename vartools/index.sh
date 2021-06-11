#!/usr/bin/bash
ref=$1
tool_name=$2
prefix=$3
soft=./softwares.config
tool=`grep $tool_name"=" $soft |cut -d '=' -f 2`
new_file=$prefix.fa
if [[ ! -e "$new_file" ]];then
    awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\n":toupper($0) }' $ref > $new_file
fi
####statistics####
genome_size=`awk '$0 !~ /^>/{print length}' $new_file |awk '{s+=$1}END{print s}'`
####index####
if [[ $tool_name = 'bwa' && $genome_size -lt 4000000000 ]] ; then  ## <4G small genome
    $tool index $new_file -p $prefix
elif [[ $tool_name = 'bwa' && $genome_size -ge 4000000000 ]] ; then ## >4G large genome
     $tool index $new_file -a bwtsw -p $prefix
#if [[ $tool_name = 'bwa2' ]]; then   ## Need too many memory, avoid to use it
#    $tool index $new_file -p $prefix
elif [[ $tool_name = 'samtools' ]] ;then
    $tool faidx $new_file
elif [[ $tool_name = 'hisat2' ]];then
    $tool"-build" -p 8 $new_file $prefix
elif [[ $tool_name = 'gatk4' ]];then
    $tool --java-options "-Xmx8G" CreateSequenceDictionary -R $new_file -O $prefix.dict
elif [[ $tool_name = 'fato2bit' ]];then
    $tool $new_file $prefix.fa.2bit

else
    echo "Error,no such tool in softwares.config or input wrong name, please add software in software.config or use correct name"
fi
