# Variation_pipelinev1.0
The common pipeline of dectecting variations:SNP/Indel/SV/CNV that combined by samtools,GATK, Breakdancer,sniffles...
```
usage: run_variation.py [-h] [-i INPUT] [-o OUTPUT] [-r REFERENCE] [-g GFF3] [-bv {hg19,hg38}] [-sp {ngs,tgs}] [-st {PE,SE}]
                        [-mt {BWA,Minimap2,NGMLR}] [-mp MAPTOOL_PARAMETERS] [-sg {WGS,WES}]
                        [-cp {samtools,gatk4,samtools+gatk4,breakdancer,crest,cnvnator,control-freec,sniffles}]
                        [-mode {SNP_INDEL,SNP_INDEL_Somatic,CNV,SV,SV_Somatic,CNV_Somatic}] [-bqsr BQSR] [-vqsr VQSR] [-vc {single,join}]
                        [-tar TARGET] [-con CONTROL] [-pon PANEL_OF_NORMALS] [-np NORMAL_SAMPLES_FOR_PON] [-gm GERMLINE]
                        [-af AF_OF_ALLELES_NOT_IN_RESOURCE] [-sniffles_p SNIFFLES_PARAMETERS]

Variation pipline v1.0

General options:
  -h, --help            show the help and exit
  -i INPUT, --input INPUT
                        The input directory of clean reads
  -o OUTPUT, --output OUTPUT
                        The output directory
  -r REFERENCE, --reference REFERENCE
                        The fasta file of reference
  -g GFF3, --gff3 GFF3  The gff3 file of reference for variation sites annotation
  -bv {hg19,hg38}, --build_version {hg19,hg38}
                        Human genome build version,if used,do not set -r and -g

Mapping options:
  -sp {ngs,tgs}, --seq_platform {ngs,tgs}
                        Reads sequencing platform, ngs for illumina short reads; tgs for Pacbio or ONT long reads
  -st {PE,SE}, --seq_type {PE,SE}
                        Sequencing type, paired-end/mate-pair, or single end
  -mt {BWA,Minimap2,NGMLR}, --maptool {BWA,Minimap2,NGMLR}
                        Choose an alignment tool,illumina:BWA;Pacbio/ONT:Minimap2/NGMLR
  -mp MAPTOOL_PARAMETERS, --maptool_parameters MAPTOOL_PARAMETERS
                        Set parameters for alignment tools

Detection variation (SNP/INDEL/SV/CNV) options:
  -sg {WGS,WES}, --strategy {WGS,WES}
                        WGS or WES
  -cp {samtools,gatk4,samtools+gatk4,breakdancer,crest,cnvnator,control-freec,sniffles}, --callpipe {samtools,gatk4,samtools+gatk4,breakdancer,crest,
cnvnator,control-freec,sniffles}
                        Choose a detection pipeline for SNP/INDEL/SV/CNV
  -mode {SNP_INDEL,SNP_INDEL_Somatic,CNV,SV,SV_Somatic,CNV_Somatic}, --mode {SNP_INDEL,SNP_INDEL_Somatic,CNV,SV,SV_Somatic,CNV_Somatic}
                        Mutation types
  -bqsr BQSR, --BQSR BQSR
                        The directory that only contain known sites for BQSR parameter(GATK4 only)
  -vqsr VQSR, --VQSR VQSR
                        Choose a method for mutations quality control(GATK4 only), if set this parameters, please provide a directory that only
                        contain resource files; if not, hard_filtering will be used
  -vc {single,join}, --variation_calling {single,join}
                        Calling a group of samples together(join calling) or Variant calling with a single sample only(single sample calling),
                        default is single
  -tar TARGET, --target TARGET
                        The target(tumor/disease) sample name
  -con CONTROL, --control CONTROL
                        The control(normal) sample name
  -pon PANEL_OF_NORMALS, --panel_of_normals PANEL_OF_NORMALS
                        The PON(panel of normal) vcf file,if not provided, use the --normal_samples_for_pon paramters to create the normal panel
  -np NORMAL_SAMPLES_FOR_PON, --normal_samples_for_pon NORMAL_SAMPLES_FOR_PON
                        The sample name of normal samples to create the PON,e.g. sample1,sample2,sample3...
  -gm GERMLINE, --germline GERMLINE
                        The population germline resource
  -af AF_OF_ALLELES_NOT_IN_RESOURCE, --af_of_alleles_not_in_resource AF_OF_ALLELES_NOT_IN_RESOURCE
                        The GATK4 af-of-alleles-not-in-resource parameter in Mutect2, The default of 0.001 is appropriate for human sample analyses
                        without any population resource.
  -sniffles_p SNIFFLES_PARAMETERS, --sniffles_parameters SNIFFLES_PARAMETERS
                        The sniffles parameters for SV calling.

EXAMPLES: 
python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -bv hg38 -sp ngs -mt BWA -cp samtools -mode SNP_INDEL
python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -bv hg19 -sp ngs -mt BWA -cp gatk4 -bqsr /root/my_data/known_sites/ -vqsr /root/my_data/resources/ -tar target_name -con control_name -mode SNP_INDEL
python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -cp gatk4 -tar Tu4 -con Nm_35 -mode SNP_INDEL
python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -cp control-freec -mode CNV -sg WES 
python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -bv hg38 -sp ngs -mt BWA -cp gatk4 -mode SNP_INDEL_Somatic -tar SRR1553607.1 -con SRR1553608
python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -bv hg38 -sp ngs -mt BWA -cp gatk4 -mode SNP_INDEL_Somatic -tar SRR1553607.1 -con SRR1553608 -np sample1,sample2,samples3 
python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -bv hg38 -sp ngs -mt BWA -cp gatk4 -mode SNP_INDEL_Somatic -tar SRR1553607.1 -con SRR1553608 -pon /root/my_data/genomicsdb/g38/hg38_pon.vcf -gm /root/my_data/germlines/hg19/af-only-gnomad.vcf.gz -af 0.001
python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g /root/my_data/ref/hg19.gff3 -sp tgs -mt Minimap2 -cp gatk4 -mode SNP_INDEL
python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g /root/my_data/ref/hg19.gff3
-sp tgs -mt NGMLR -cp sniffles -mode SV

```
