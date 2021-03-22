import os
from tools import getmyconfig

samtools = getmyconfig.getConfig('Variation', 'samtools')
bcftools = getmyconfig.getConfig('Variation', 'bcftools')
gatk = getmyconfig.getConfig('Variation', 'GATK')
breakdancer = getmyconfig.getConfig('Variation', 'BreakDancer')
cnvnator = getmyconfig.getConfig('Variation', 'CNVnator')
sniffles = getmyconfig.getConfig('Variation', 'Sniffles')
annovar = getmyconfig.getConfig('Variation', 'ANNOVAR')

def ngs_snp_indel(ref, sample, tool, parameters, combine='N'):  ## Samtools and GATK pipelines
    outfile = ''
    if combine == "N" and tool == 'samtools':
        outfile = """{samtools} {parameters} {ref} {sample}.rmdup.bam | {bcftools} call -m -f GQ | gzip > {sample}.all.vcf.gz
        python {getVarformAll} -i {sample}.all.vcf.gz -o {sample}.raw.vcf
        perl {vcfutils} varFilter -Q 20 -d 4 -D 1000 {sample}.raw.vcf > {sample}.filter.vcf
        python {splitSNPindelVCF} -i {sample}.filter.vcf -p {sample}.filted
        """.format(samtools=tool,parameters=parameters,ref=ref,sample=sample,bcftools=bcftools,getVarformAll=getVarformAll,
                   vcfutils=vcfutils,splitSNPindelVCF=splitSNPindelVCF)
    elif combine == 'N' and tool == 'gatk':
        outfile = """jar {}"""
    elif combine == 'Y' and tool == 'samtools+gatk':


def ngs_sv(file, tool, parameters):

def ngs_cnv(file, tool, parameters):
