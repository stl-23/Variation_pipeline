# Variation_pipelinev1.0
The common pipeline of dectecting variations:SNP/Indel/SV/CNV that combined by samtools,GATK, Breakdancer,sniffles...
# Dependencies download and usage
1. Pull our docker image from Docker Hub 
```
docker pull stl23/variation:v1.0
```
2. Download reference genome data (hg19 or hg38) 

Note:ANNOVAR is used for annotation,please provide the path of ANNOVAR software.
```
wget -c https://github.com/stl-23/Variation_pipelinev1.0/install_genome_data.sh

cd “YOUR_DATABASE_DIR” && sh ./install_genome_data.sh hg38 /tool/annovar
```
3 Start docker image
docker run -itd --name variationv1.0 -v "YOUR_DATABASE_DIR":/scripts/database/genomicsdb/ -v “YOUR_WORK_DIR”:/work -v "YOUR_INPUT_DIR":/input -v "YOUR_OUTPUT_DIR":/output variation:v1.0 /bin/bash
```
4. Run the pipeline
```
docker run -v "YOUR_INPUT_DIR":/input\
-v "YOUR_OUTPUT_DIR":/output \
-v “YOUR_WORK_DIR”:/work \
variation:v1.0 bash -c 'cd /work/ && python3 /script/variation/run_variation.py \
-i /input \
-o /output \
-bv hg38 \
-sp ngs \
-sg WGS \
-mt BWA \
-cp gatk4 \
-mode SNP_INDEL
```
# Parameters
```
usage: run_variation.py [-h] [-i INPUT] [-o OUTPUT] [-r REFERENCE] [-g GFF3]
                        [-bv {hg19,hg38}] [-mx MAX_MEMORY] [-mc MAX_THREAD]
                        [-sp {ngs,tgs}] [-st {PE,SE}]
                        [-mt {BWA,Minimap2,NGMLR}] [-mp MAPTOOL_PARAMETERS]
                        [-sg {WGS,WES}]
                        [-cp {samtools,gatk4,samtools+gatk4,breakdancer,cnvnator,control-freec,sniffles}]
                        [-mode {SNP_INDEL,SNP_INDEL_Somatic,CNV,SV,CNV_Somatic}]
                        [--bcftools_filter BCFTOOLS_FILTER] [-bqsr [BQSR]]
                        [-vqsr [VQSR]] [-vc {single,join}] [-tar TARGET]
                        [-con CONTROL] [-pon PANEL_OF_NORMALS]
                        [-np NORMAL_SAMPLES_FOR_PON] [-gm GERMLINE]
                        [-af AF_OF_ALLELES_NOT_IN_RESOURCE]
                        [--interval INTERVAL]
                        [--sniffles_parameters SNIFFLES_PARAMETERS]

Variation pipeline v1.0

General options:
  -h, --help            show the help and exit
  -i INPUT, --input INPUT
                        The input directory of clean reads
  -o OUTPUT, --output OUTPUT
                        The output directory
  -r REFERENCE, --reference REFERENCE
                        The fasta file of reference
  -g GFF3, --gff3 GFF3  The gff3 file of reference for variation sites
                        annotation
  -bv {hg19,hg38}, --build_version {hg19,hg38}
                        Human genome build version,if used,do not set -r and
                        -g
  -mx MAX_MEMORY, --max_memory MAX_MEMORY
                        The maximum JAVA memory in the pipeline
  -mc MAX_THREAD, --max_thread MAX_THREAD
                        The maximum Java theads for GATK GenomicsDBImport

Mapping options:
  -sp {ngs,tgs}, --seq_platform {ngs,tgs}
                        Reads sequencing platform, ngs for illumina short
                        reads; tgs for Pacbio or ONT long reads
  -st {PE,SE}, --seq_type {PE,SE}
                        Sequencing type, paired-end/mate-pair, or single end
  -mt {BWA,Minimap2,NGMLR}, --maptool {BWA,Minimap2,NGMLR}
                        Choose an alignment
                        tool,illumina:BWA;Pacbio/ONT:Minimap2/NGMLR
  -mp MAPTOOL_PARAMETERS, --maptool_parameters MAPTOOL_PARAMETERS
                        Set parameters for alignment tools

Detection variation (SNP/INDEL/SV/CNV) options:
  -sg {WGS,WES}, --strategy {WGS,WES}
                        WGS or WES
  -cp {samtools,gatk4,samtools+gatk4,breakdancer,cnvnator,control-freec,sniffles}, --callpipe {samtools,gatk4,samtools+gatk4,breakdancer,cnvnator,control-freec,sniffles}
                        Choose a detection pipeline for SNP/INDEL/SV/CNV
  -mode {SNP_INDEL,SNP_INDEL_Somatic,CNV,SV,CNV_Somatic}, --mode {SNP_INDEL,SNP_INDEL_Somatic,CNV,SV,CNV_Somatic}
                        Mutation types
  --bcftools_filter BCFTOOLS_FILTER
                        Hard filter paramters for SNP/Indel detection pipeline
                        by bcftools. See
                        http://samtools.github.io/bcftools/howtos/variant-
                        calling.html
  -bqsr [BQSR], --BQSR [BQSR]
                        The directory that only contain known sites for BQSR
                        parameter(GATK4 only)
  -vqsr [VQSR], --VQSR [VQSR]
                        Choose a method for mutations quality control(GATK4
                        only), if set this parameters, please provide a
                        directory that only contain resource files; if not,
                        hard_filtering will be used
  -vc {single,join}, --variation_calling {single,join}
                        Calling a group of samples together(join calling) or
                        Variant calling with a single sample only(single
                        sample calling), default is single

Somatic mutation options:
  -tar TARGET, --target TARGET
                        The target(case/tumor/disease) sample name
  -con CONTROL, --control CONTROL
                        The control(normal) sample name
  -pon PANEL_OF_NORMALS, --panel_of_normals PANEL_OF_NORMALS
                        The PON(panel of normal) vcf file,if not provided, use
                        the --normal_samples_for_pon paramters to create the
                        normal panel
  -np NORMAL_SAMPLES_FOR_PON, --normal_samples_for_pon NORMAL_SAMPLES_FOR_PON
                        The sample name of normal samples to create the
                        PON,e.g. sample1,sample2,sample3...
  -gm GERMLINE, --germline GERMLINE
                        The population germline resource
  -af AF_OF_ALLELES_NOT_IN_RESOURCE, --af_of_alleles_not_in_resource AF_OF_ALLELES_NOT_IN_RESOURCE
                        The GATK4 af-of-alleles-not-in-resource parameter in
                        Mutect2, The default of 0.001 is appropriate for human
                        sample analyses without any population resource.
  --interval INTERVAL   The interval file

TGS(Third generation sequencing) SV options:
  --sniffles_parameters SNIFFLES_PARAMETERS
                        The sniffles parameters for SV calling.
```
# EXAMPLES:
```
## samtools single calling
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -sp ngs -sg WGS -mt BWA -cp samtools -mode SNP_INDEL
## gatk4 single calling
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -sp ngs -sg WGS -mt BWA -cp gatk4 -mode SNP_INDEL
## join calling
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -st PE -sp ngs -sg WGS -mt BWA -cp samtools -vc join -mode SNP_INDEL
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -st PE -sp ngs -sg WGS -mt BWA -cp gatk4 -vc join -mode SNP_INDEL
## bqsr & vqsr
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -st PE -sp ngs -sg WGS -mt BWA -cp gatk4 -vc join -bqsr -vqsr -mode SNP_INDEL
## samtools+gatk4
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -st PE -sp ngs -sg WGS -mt BWA -cp samtools+gatk4 -vc join -bqsr -vqsr -mode SNP_INDEL
## somatic SNP/Indel
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -st PE -sp ngs -sg WGS -mt BWA -cp gatk4 -tar H9-AB-p113_S2_L001 -con H9-AB-p113_S10_L001 -mode SNP_INDEL_Somatic
## somatic SNP/Indel create PON
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -st PE -sp ngs -sg WGS -mt BWA -cp gatk4 -tar H9-AB-p113_S2_L001 -con H9-AB-p113_S10_L001 -np H9-AB-p113_S13_L001,H9-AB-p113_S5_L001,H9-AB-p116_S4_L001 -mode SNP_INDEL_Somatic
## SV
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -st PE -sp ngs -mt BWA -cp breakdancer -mode SV'
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -st PE -sp ngs -mt BWA -cp crest -mode SV
## CNV
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -st PE -sg WGS -sp ngs -mt BWA -cp cnvnator -mode CNV
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -st PE -sg WGS -sp ngs -mt BWA -cp control-freec -mode CNV
## somatic CNV 
python3 /scripts/variation/run_variation.py -i /qc/results/Illumina/ -o /variation/results/ -bv hg38 -st PE -sg WGS -sp ngs -mt BWA -cp control-freec -tar H9-AB-p113_S2_L001 -con H9-AB-p113_S10_L001 -mode CNV_Somatic
## TGS SV 
python3 /scripts/variation/run_variation.py -i /qc/results/pacbio/DNA/ -o /variation/results/ -bv hg38 -sp tgs -mt Minimap2 -mp "--cs --MD -ax map-pb" -cp sniffles -mode SV
python3 /scripts/variation/run_variation.py -i /qc/results/pacbio/DNA/ -o /variation/results/ -bv hg38 -sp tgs -mt NGMLR -mp "-x pacbio" -cp sniffles -mode SV
```
