#!/usr/bin/env bash
if [[ "$#" -lt 1 ]]; then
  echo
  echo "This script downloads and installs data for genomes."
  echo "Genome data files will be installed on local path ./"
  echo
  echo "Supported genomes: hg19, hg38"
  echo
  echo "Usage: ./install_genome_data.sh [GENOME]"
  echo "  Example: ./install_genome_data.sh hg19"
  echo
  exit 0
fi

GENOME=$1
Annovar_path=$2

if [[ $GENOME == "hg19" ]]; then
  #ref_fa="ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz"
  ref_fa="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta"
  gtf="ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz"
  dict="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict"
  amb="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.amb"
  ann="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.ann"
  bwt="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.bwt"
  fai="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai"
  pac="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.pac"
  sa="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.sa"
  ## resource files
  hapmap="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/hapmap_3.3.b37.vcf.gz"
  hapmap_index="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/hapmap_3.3.b37.vcf.gz.tbi"
  omini="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/1000G_omni2.5.b37.vcf.gz"
  omini_index="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/1000G_omni2.5.b37.vcf.gz.tbi"
  G1000="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz"
  G1000_index="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz.tbi"
  dbsnp="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz"
  dbsnp_index="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz.tbi"
  mills="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
  mills_index="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Mills_and_1000G_gold_standard.indels.b37.vcf.gz.tbi"
  axiomPoly="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz"
  axiomPoly_index="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz.tbi"
  interval="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/wgs_calling_regions.v1.interval_list"
  ## bqsr resource
  known_indel="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.known_indels.vcf"
  known_indel_idx="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.known_indels.vcf.idx"

elif [[ $GENOME == "hg38" ]]; then
  #ref_fa="ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz"
  ref_fa="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
  gtf="ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz"
  dict="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
  amb="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.amb"
  ann="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.ann"
  bwt="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.bwt"
  fai="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
  pac="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.pac"
  sa="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.sa"

  hapmap="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz"
  hapmap_index="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi"
  omini="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz"
  omini_index="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi"
  G1000="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  G1000_index="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
  dbsnp="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
  dbsnp_index="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
  mills="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  mills_index="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
  axiomPoly="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
  axiomPoly_index="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"

  interval="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list"
  known_indel="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
  known_indel_idx="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
fi

if [[ ${ref_fa} == "" ]]; then
  echo "Error: unsupported genome $GENOME"
  exit 1
fi
echo

mkdir ./vqsr_resource ./bqsr_resource ./pon_and_germline ./annovar ./split_fa/
## download
for i in $ref_fa $gtf $dict $amb $ann $bwt $fai $pac $sa
do
  {
    wget -N -c $i
  }&
done

for i in $hapmap $hapmap_index $omini $omini_index $G1000 $G1000_index
do
  {
    wget -N -c $i
  }&
done

for i in $dbsnp $dbsnp_index $mills $mills_index $axiomPoly $axiomPoly_index $interval
do
  {
    wget -N -c $i
  }&
done

for i in $known_indel $known_indeld_idx
do
  {
    wget -N -c $i
  }&
## rename
rename Homo_sapiens_assembly hg Homo_sapiens_assembly*
rename fasta fa *fasta*
if [[ $GENOME == "hg19" ]]; then
  ## resource
  mv *.vcf.gz *.tbi ./vqsr_resource
  mv Homo_sapiens_assembly19.known_indels.vcf* ./bqsr_resource
  cd ./bqsr_resource && ln -s ../vqsr_resource/dbsnp_138.b37.vcf.gz* ./ && ln -s ../vqsr_resource/Mills_and_1000G_gold_standard.indels.b37.vcf.gz* ./
  mv wgs_calling_regions.v1.interval_list hg19.interval_list
  cd ../
  ## fasta file
  rename Homo_sapiens_assembly19 hg19 Homo_sapiens_assembly19*
  rename fasta fa *fasta
  gunzip hg19.refGene.gtf.gz && mv hg19.refGene.gtf hg19.gtf
  ## split fasta file
  cd ./split_fa
  for i in `grep '>' ../hg19.fa |sed 's/>//g'`
  do
    echo $i > $i.txt && awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' $i.txt ../hg19.fa > $i.fa && rm $i.txt
  done
  cd ../
  ## annovar
  cd ./annovar
  perl $Annovar_path/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/

elif [[ $GENOME == "hg38" ]]; then
    ## resource
  mv *.vcf.gz *.tbi ./vqsr_resource
  mv Homo_sapiens_assembly38.known_indels.vcf.gz* ./bqsr_resource
  cd ./bqsr_resource && ln -s ../vqsr_resource/Homo_sapiens_assembly38.dbsnp138.vcf.gz* ./ && ln -s ../vqsr_resource/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz* ./
  mv wgs_calling_regions.hg38.interval_list hg38.interval_list
  cd ../
  ## fasta file
  rename Homo_sapiens_assembly38 hg38 Homo_sapiens_assembly38*
  rename fasta fa *fasta
  gunzip hg38.refGene.gtf.gz && mv hg38.refGene.gtf hg38.gtf
  ## split fasta file
  cd ./split_fa
  for i in `grep '>' ../hg38.fa |sed 's/>//g'`
  do
    echo $i > $i.txt && awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' $i.txt ../hg38.fa > $i.fa && rm $i.txt
  done
  cd ../
  ## annovar
  cd ./annovar
  perl $Annovar_path/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
