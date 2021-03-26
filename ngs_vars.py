import os,sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath('./'))))
from vartools import getmyconfig,make_freec_config

samtools = getmyconfig.getConfig('Variation', 'samtools')
bcftools = getmyconfig.getConfig('Variation', 'bcftools')
vcfutils = getmyconfig.getConfig('Variation','vcfutils')
gatk4 = getmyconfig.getConfig('Variation', 'gatk4')
breakdancer = getmyconfig.getConfig('Variation', 'BreakDancer')
bam2cfg = getmyconfig.getConfig('Variation','bam2cfg')
crest = getmyconfig.getConfig('Variation','Crest')
extractSClip = getmyconfig.getConfig('Variation','extractSClip')
cnvnator = getmyconfig.getConfig('Variation', 'CNVnator')
cnvnator2VCF = getmyconfig.getConfig('Variation','cnvnator2VCF')
control_freec = getmyconfig.getConfig('Variation','control_freec')
freec_WGS_config = getmyconfig.getConfig('Variation','freec_WGS_config')
freec_WES_config = getmyconfig.getConfig('Variation','freec_WES_config')


## known_site='--known-sites /path/to/ref1.vcf --known-sites /path/to/ref2.vcf ....'
def snp_indel_samtools(ref, sample):
    outfile = ''
    samtools_p = 'mpileup -C 50 -m 2 -F 0.002 -d 1000 -u -f'  ### Hard filtering
    bcftools_p = 'call -m -f GQ'
    vcfutils_p = 'varFilter -Q 20 -d 4 -D 1000'
    outfile = """{samtools} {samtools_p} {ref} {sample}.rmdup.bam | {bcftools} {bcftools_p} | gzip > {sample}.all.vcf.gz
        python {getVarformAll} -i {sample}.all.vcf.gz -o {sample}.raw.vcf
        perl {vcfutils} {vcfutils_p} {sample}.raw.vcf > {sample}.filter.vcf
        python {splitSNPindelVCF} -i {sample}.filter.vcf -p {sample}.samtools.raw
        """.format(samtools=tool,samtools_p=samtools_p,ref=ref,sample=sample,bcftools=bcftools,getVarformAll=getVarformAll,
                   vcfutils=vcfutils,splitSNPindelVCF=splitSNPindelVCF,bcftools_p=bcftools_p,vcfutils_p=vcfutils_p)
    return outfile

def snp_indel_gatk(ref, sample, BQSR='T', gvcf='T',known_site):
    outfile = ''
    if BQSR == "T":  #### BQSR (Recalibration Base Quality Score)
        if gvcf == "T":
            #if filter == "VQSR":
            outfile = """{gatk4} BaseRecalibrator -R {ref} -I {sample}.rmdup.bam {known_site} -O {sample}.recal.table
                    {gatk4} ApplyBQSR --bqsr-recal-file {sample}.recal.table -R {ref} -I {sample}.rmdup.bam -O {sample}.BQSR.bam
                    {gatk4} HaplotypeCaller -R {ref} -I {sample}.BQSR.bam -ERC GVCF -O {sample}.g.vcf
                    """.format(gatk4=gatk4, ref=ref, sample=sample, known_site=known_site)
            #elif filter == "Hard_filter":

        elif gvcf == "F":
           # if filter == "VQSR":
            outfile = """{gatk4} BaseRecalibrator -R {ref} -I {sample}.rmdup.bam {known_site} -O {sample}.recal.table
                       {gatk4} ApplyBQSR --bqsr-recal-file {sample}.recal.table -R {ref} -I {sample}.rmdup.bam -O {sample}.BQSR.bam
                       {gatk4} HaplotypeCaller -R {ref} -I {sample}.BQSR.bam -O {sample}.vcf
                        """.format(gatk4=gatk4, ref=ref, sample=sample, known_site=known_site)
    elif BQSR == "F":
        if gcvf == "T":
            outfile = """{gatk4} HaplotypeCaller -R {ref} -I {sample}.BQSR.bam -ERC GVCF -O {sample}.g.vcf""".format(gatk4=gatk4, ref=ref, sample=sample)
        elif gcvf == "F":
            outfile = """{gatk4} HaplotypeCaller -R {ref} -I {sample}.BQSR.bam -O {sample}.vcf""".format(gatk4=gatk4, ref=ref, sample=sample)
    return outfile
'''
def ngs_snp_indel(ref, sample, tool, combine='N',BQSR,gcvf,known_site):  ## Samtools and GATK pipelines
    out1 = ''
    out2 = ''
    out3 = ''
    if combine == "N" and tool == 'samtools':
        out1 = snp_indel_samtools(ref, sample)
    elif combine == 'N' and tool == 'gatk':
        out2 = snp_indel_gatk(ref, sample, BQSR, gcvf,known_site)
    elif combine == 'Y' and tool == 'samtools+gatk':
        out1 = snp_indel_samtools(ref, sample)
        out2 = snp_indel_gatk(ref, sample, BQSR, gcvf,known_site)
        out3 = """{gatk4} SelectVariants --variant {sample}.gatk.raw.vcf --concordance {sample}.samtools.raw.vcf -o {sample}.concordance.raw.vcf""".format(
           gatk4=gatk4,sample=sample
        )
    return out1,out2,out3
'''
def ngs_sv(sample1,sample2,ref,tool='breakdancer',rm_germline="false"): ## sample1:disease/somatic; sample2:control
    outfile = ''
    breakdancer_p = '-q 20 -d'
    #crest_p = ''
    if tool == 'breakdancer':  ## Only use paired end reads
        if rm_germline == "true": ## somatic SV
            pass
        elif rm_germline == "false": ## germline/common SV
            outfile = """perl {bam2cfg} q 20 -c 4 -g -h {sample1}.rmdup.bam {sample1}.bam.cfg
                {breakdancer} {breakdancer_p} {sample1} {sample1}.bam.cfg > {sample1}.raw.ctx""".format(
            bam2cfg=bam2cfg,sample1=sample1,breakdancer=breakdancer,breakdancer_p=breakdancer_p
        )
    elif tool == 'crest':      ## Can use both paired end and singe end reads
        if rm_germline == "true":  ## somatic SV
            outfile = """perl {extractSClip} -i {sample1}.rmdup.bam --ref_genome {ref} -p {sample1}
            perl {crest} -f {sample1}.cover -d {sample1}.rmdup.bam -g {sample2}.rmdup.bam --ref_genome {ref} -t {ref}.2bit -p {sample1}
            """.format(
                extractSClip=extractSClip,sample1=sample1,sample2=sample2,
                crest=crest,ref=ref
            )
        elif rm_germline == "false": ## germline/common SV
            outfile = f"""perl {extractSClip} -i {sample1}.rmdup.bam --ref_genome {ref} -p {sample1}
             perl {crest} -f {sample1}.cover -d {sample1}.rmdup.bam --ref_genome {ref} -t {ref}.2bit -p {sample1}""".format(
                extractSClip=extractSClip, sample1=sample1,crest=crest,ref=ref
            )
    return outfile
def ngs_cnv(sample1, sample2, ref, tool='control-freec',species='human',omic_type="WGS",rm_germline="false",outdir):
    outfile = ''
    if tool == 'cnvnator':
        bin_size = '1000'
        ref_dir = os.path.dirname(ref)
        outfile = """{cnvnator} -root {sample1}.root -tree {sample1}.rmdup.bam 
        {cnvnator} -root {sample1}.root -his {bin_size} -d {ref_dir}
        {cnvnator} -root {sample1}.root -stat {bin_size}
        {cnvnator} -root {sample1}.root -partition {bin_size}
        {cnvnator} -root {sample1}.root -call {bin_size} > {sample1}.cnv.call.txt
        perl {cnvnator2VCF} {sample1}.cnv.call.txt > {sample1}.cnv.vcf
        """.format(cnvnator=cnvnator,sample1=sample1,ref_dir=ref_dir,cnvnator2VCF=cnvnator2VCF,bin_size=bin_size)
    elif tool == 'control-freec':
        if species == 'human':
            if omic_type == "WGS":
                if rm_germline == "true":
                    sample_data = sample1+'rmdup.bam'
                    control_data = sample2+'rmdup.bam'
                    config_WGS_add_control.txt=make_freec_config.modify(sample_data,control_data,'human','WGS','Y',outdir)
                    outfile = """{control_freec} -conf {config_WGS_add_control}""".format(control_freec=control_freec,
                                                                    config_WGS_add_control=config_WGS_add_control.txt)
                elif rm_germline == "false":
                    sample_data = sample1+'rmdup.bam'
                    config_WGS_no_control.txt=make_freec_config.modify(sample_data,'','human','WGS','N',outdir)
                    outfile = """{control_freec} -conf {config_WGS_no_control}""".format(control_freec=control_freec,
                                                                    config_WGS_no_control=config_WGS_no_control.txt)
            elif omic_type == "WES":
                if rm_germline == "true":
                    sample_data = sample1 + 'rmdup.bam'
                    control_data = sample2 + 'rmdup.bam'
                    config_WES_add_control.txt = make_freec_config.modify(sample_data,control_data,'human','WES','Y',outdir)
                    outfile = """{control_freec} -conf {config_WES_add_control}""".format(control_freec=control_freec,
                                                                    config_WES_add_control=config_WES_add_control.txt)
                elif rm_germline == "false":
                    sample_data = sample1 + 'rmdup.bam'
                    config_WES_add_control.txt = make_freec_config.modify(sample_data, '', 'human', 'WES', 'N', outdir)
                    outfile = """{control_freec} -conf {config_WES_add_control}""".format(control_freec=control_freec,
                                                           config_WES_add_control=config_WES_add_control.txt)
        elif species == 'non-human':
            pass

    elif tool == 'ExomeCNV': ## WES
        pass