
from vartools import getmyconfig

samtools = getmyconfig.getConfig('Variation', 'samtools')
sniffles = getmyconfig.getConfig('Variation', 'Sniffles')

def tgs_snp_indel():
    pass

def tgs_sv(sample,tool='sniffles',sniffle_p='-s 1 -d 600 --genotype --cluster --ccs_reads'):
    outfile = ''
    if tool == 'sniffles':
        outfile = """{sniffles} {sniffle_p} -m {sample}.sorted.bam -v {sample}.sv.vcf""".format(
            sniffles=sniffles,sniffle_p=sniffle_p,sample=sample
        )
    return outfile