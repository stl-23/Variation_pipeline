
from vartools import getmyconfig

gatk4 = getmyconfig.getConfig('Variation', 'gatk4')
samtools = getmyconfig.getConfig('Variation', 'samtools')
sniffles = getmyconfig.getConfig('Variation', 'Sniffles')

def tgs_snp_indel(ref,input,sample):
    ### gatk4 pipeline ###
    outfile = """{gatk4} HaplotypeCaller -R {ref} -I {input}.sorted.bam \\
 --pcr_indel_model AGGRESSIVE \\
 --annotation-group AS_StandardAnnotation \\
 --minimum-mapping-quality 60 \\
-O {sample}.vcf

{gatk4} SelectVariants  -R {ref} -V {sample}.vcf --select-type-to-include SNP -O {sample}.raw.snp.vcf
{gatk4} SelectVariants  -R {ref} -V {sample}.vcf --select-type-to-include INDEL -O {sample}.raw.indel.vcf

{gatk4} VariantFiltration -R {ref} -V {sample}.raw.snp.vcf --filter-expression \\

""".format(gatk4=gatk4, ref=ref, input=input,sample=sample)
    return outfile

def tgs_sv(sample,tool='sniffles',sniffle_p='-s 1 -d 600 --genotype --cluster --ccs_reads'):
    outfile = ''
    if tool == 'sniffles':
        outfile = """{sniffles} {sniffle_p} -m {sample}.sorted.bam -v {sample}.sv.vcf""".format(
            sniffles=sniffles,sniffle_p=sniffle_p,sample=sample
        )
    return outfile