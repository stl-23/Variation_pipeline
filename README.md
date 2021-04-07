# Variation_pipelinev1.0
The common pipeline of dectecting variations:SNP/Indel/SV/CNV that combined by samtools,GATK, Breakdancer,sniffles...
```
usage: run_variation.py [-h] [-i INPUT] [-o OUTPUT] [-r REFERENCE] [-g GFF3] [-sp {ngs,tgs}] [-st {PE,SE}] [-mt {BWA,Minimap2,NGML}]
                        [-mp MAPTOOL_PARAMETERS] [-sg {WGS,WES}]
                        [-cp {samtools,gatk4,samtools+gakt4,breakdancer,crest,cnvnator,control-freec,sniffles}]
                        [-mode {SNP_Indel,SV,CNV,Germline,Somatic}] [-bqsr BQSR] [-vqsr VQSR] [-bv {hg19,hg38,other}] [-vc {single,join}]
                        [-s1 TARGET] [-s2 CONTROL] [-pon PANEL_OF_NORMALS] [-np NORMAL_SAMPLES_FOR_PON] [-gm GERMLINE]
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

Mapping options:
  -sp {ngs,tgs}, --seq_platform {ngs,tgs}
                        Reads sequencing platform, ngs for illumina short reads; tgs for Pacbio or ONT long reads
  -st {PE,SE}, --seq_type {PE,SE}
                        Sequencing type, paired-end/mate-pair, or single end
  -mt {BWA,Minimap2,NGML}, --maptool {BWA,Minimap2,NGML}
                        Choose an alignment tool,illumina:BWA;Pacbio/ONT:Minimap2/NGML
  -mp MAPTOOL_PARAMETERS, --maptool_parameters MAPTOOL_PARAMETERS
                        Set parameters for alignment tools

Detection variation (SNP/Indel/SV/CNV) options:
  -sg {WGS,WES}, --strategy {WGS,WES}
                        WGS or WES
  -cp {samtools,gatk4,samtools+gakt4,breakdancer,crest,cnvnator,control-freec,sniffles}, --callpipe {samtools,gatk4,samtools+gakt4,breakdancer,crest,
cnvnator,control-freec,sniffles}
                        Choose a detection pipeline for SNP/Indel/SV/CNV
  -mode {SNP_Indel,SV,CNV,Germline,Somatic}, --mode {SNP_Indel,SV,CNV,Germline,Somatic}
                        Mutation types
  -bqsr BQSR, --BQSR BQSR
                        If choose BQSR parameter(GATK4 only), you should provide a directory that only contain known sites
  -vqsr VQSR, --VQSR VQSR
                        Choose a method for mutations quality control(GATK4 only), if set this parameters, please provide a directory that only
                        contain resource files, if not hard_filtering will be used
  -bv {hg19,hg38,other}, --build_version {hg19,hg38,other}
                        Species and genome version
  -vc {single,join}, --variation_calling {single,join}
                        Calling a group of samples together(join calling) or Variant calling with a single sample only(single sample calling),
                        default is single
  -s1 TARGET, --target TARGET
                        The target(tumor/disease) sample name
  -s2 CONTROL, --control CONTROL
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
```
