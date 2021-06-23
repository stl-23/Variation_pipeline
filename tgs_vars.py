
from vartools import getmyconfig

gatk4 = getmyconfig.getConfig('Variation', 'gatk4')
samtools = getmyconfig.getConfig('Variation', 'samtools')
sniffles = getmyconfig.getConfig('Variation', 'Sniffles')

def tgs_snp_indel(ref,input,sample):
    ### gatk4 pipeline ###
    input_new = input.replace('.rmdup.bam', '')
    outfile = """{gatk4} HaplotypeCaller -R {ref} -I {input}.sorted.bam \\
 --pcr-indel-model AGGRESSIVE \\
 --annotation-group AS_StandardAnnotation \\
 --minimum-mapping-quality 60 \\
-O {sample}.vcf

{gatk4} SelectVariants  -R {ref} -V {sample}.vcf --select-type-to-include SNP -O {sample}.raw.snp.vcf
{gatk4} SelectVariants  -R {ref} -V {sample}.vcf --select-type-to-include INDEL -O {sample}.raw.indel.vcf

{gatk4} VariantFiltration -R {ref} -V {sample}.raw.snp.vcf \\
-filter "AS_QD < 2.0" --filter-name "ASQD2" \\
-O {sample}.snps.gatk.vcf

{gatk4} VariantFiltration -R {ref} -V {sample}.raw.indel.vcf \\
-filter "AS_QD < 5.0" --filter-name "ASQD5" \\
-O {sample}.indel.gatk.vcf
""".format(gatk4=gatk4, ref=ref, input=input_new,sample=sample)
    return outfile

def tgs_sv(sample,out_sample,tool='sniffles',sniffle_p='-s 1 -d 600 --genotype --cluster --ccs_reads'):
    outfile = ''
    sample_new = sample.replace('.rmdup.bam','')
    if tool == 'sniffles':
        outfile = """{sniffles} {sniffle_p} -m {sample}.sorted.bam -v {out}.sv.vcf""".format(
            sniffles=sniffles,sniffle_p=sniffle_p,sample=sample_new,out=out_sample
        )
    return outfile