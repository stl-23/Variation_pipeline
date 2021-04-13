import os,sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath('./'))))
from vartools import getmyconfig,make_freec_config

samtools = getmyconfig.getConfig('Variation', 'samtools').strip("'")
bcftools = getmyconfig.getConfig('Variation', 'bcftools').strip("'")
vcfutils = getmyconfig.getConfig('Variation','vcfutils').strip("'")
gatk4 = getmyconfig.getConfig('Variation', 'gatk4').strip("'")
breakdancer = getmyconfig.getConfig('Variation', 'BreakDancer').strip("'")
bam2cfg = getmyconfig.getConfig('Variation','bam2cfg').strip("'")
crest = getmyconfig.getConfig('Variation','Crest').strip("'")
extractSClip = getmyconfig.getConfig('Variation','extractSClip').strip("'")
cnvnator = getmyconfig.getConfig('Variation', 'CNVnator').strip("'")
cnvnator2VCF = getmyconfig.getConfig('Variation','cnvnator2VCF').strip("'")
control_freec = getmyconfig.getConfig('Variation','control_freec').strip("'")
freec_WGS_config = getmyconfig.getConfig('Variation','freec_WGS_config').strip("'")
freec_WES_config = getmyconfig.getConfig('Variation','freec_WES_config').strip("'")
getVarformAll = getmyconfig.getConfig('Variation','getVarformAll').strip("'")
splitSNPindelVCF = getmyconfig.getConfig('Variation','splitSNPindelVCF').strip("'")

## known_site='--known-sites /path/to/ref1.vcf --known-sites /path/to/ref2.vcf ....'
def snp_indel_samtools(ref, input, sample, v_valling):
    outfile = ''
    samtools_p = 'mpileup -C 50 -m 2 -F 0.002 -d 1000 -u -f'  ### Hard filtering
    bcftools_p = 'call -m -f GQ'
    vcfutils_p = 'varFilter -Q 20 -d 4 -D 1000'
    if v_valling == 'single':  ## input is a single bam file
        outfile = """{samtools} {samtools_p} {ref} {input_bam} | {bcftools} {bcftools_p} | gzip > {sample}.all.vcf.gz
python {getVarformAll} -i {sample}.all.vcf.gz -o {sample}.raw.vcf
perl {vcfutils} {vcfutils_p} {sample}.raw.vcf > {sample}.filter.vcf
python {splitSNPindelVCF} -i {sample}.filter.vcf -p {sample}.samtools.raw
        """.format(samtools=samtools,samtools_p=samtools_p,ref=ref,sample=sample,input_bam=input,bcftools=bcftools,
                   getVarformAll=getVarformAll,
                   vcfutils=vcfutils,splitSNPindelVCF=splitSNPindelVCF,bcftools_p=bcftools_p,vcfutils_p=vcfutils_p)
    elif v_valling == 'join': ## input is a list of bam files
        input_bams = ' '.join(input)
        outfile = """{samtools} {samtools_p} {ref} {input_bams} | {bcftools} {bcftools_p} | gzip > {sample}.all.vcf.gz
python {getVarformAll} -i {sample}.all.vcf.gz -o {sample}.raw.vcf
perl {vcfutils} {vcfutils_p} {sample}.raw.vcf > {sample}.filter.vcf
python {splitSNPindelVCF} -i {sample}.filter.vcf -p {sample}.samtools.raw""".format(
            samtools=samtools, samtools_p=samtools_p, ref=ref, sample=sample, input_bams=input_bams, bcftools=bcftools,
            getVarformAll=getVarformAll,
            vcfutils=vcfutils, splitSNPindelVCF=splitSNPindelVCF, bcftools_p=bcftools_p, vcfutils_p=vcfutils_p
        )
    return outfile

def snp_indel_gatk(ref, input, sample, gvcf, bqsr_dir):
    lst = os.listdir(os.path.abspath(bqsr_dir) + '/')
    known_site = ["--known_site "+vcf for vcf in lst if vcf.endswith('.vcf')]
    outfile = ''
    if bqsr_dir:  #### BQSR (Recalibration Base Quality Score)
        if gvcf == "T":
            outfile = """{gatk4} BaseRecalibrator -R {ref} -I {input}.rmdup.bam {known_site} -O {sample}.recal.table
                    {gatk4} ApplyBQSR --bqsr-recal-file {sample}.recal.table -R {ref} -I {sample}.rmdup.bam -O {sample}.BQSR.bam
                    {gatk4} HaplotypeCaller -R {ref} -I {sample}.BQSR.bam -ERC GVCF -O {sample}.g.vcf
                    """.format(gatk4=gatk4, ref=ref, input=input, sample=sample, known_site=known_site)

        elif gvcf == "F":
            outfile = """{gatk4} BaseRecalibrator -R {ref} -I {input}.rmdup.bam {known_site} -O {sample}.recal.table
                       {gatk4} ApplyBQSR --bqsr-recal-file {sample}.recal.table -R {ref} -I {sample}.rmdup.bam -O {sample}.BQSR.bam
                       {gatk4} HaplotypeCaller -R {ref} -I {sample}.BQSR.bam -O {sample}.vcf
                        """.format(gatk4=gatk4, ref=ref, input=input, sample=sample, known_site=known_site)
    elif not bqsr_dir:
        if gvcf == "T":
            outfile = """{gatk4} HaplotypeCaller -R {ref} -I {input}.rmdup.bam -ERC GVCF -O {sample}.g.vcf""".format(gatk4=gatk4, ref=ref, input=input, sample=sample)
        elif gvcf == "F":
            outfile = """{gatk4} HaplotypeCaller -R {ref} -I {input}.rmdup.bam -O {sample}.vcf""".format(gatk4=gatk4, ref=ref, input=input, sample=sample)
    return outfile

def samtool_gatk_combine(sample):  ## Samtools and GATK pipelines
    cmd = """{gatk4} SelectVariants --variant {sample}.snp.indel.genotype.selected.gatk.vcf --concordance {sample}.samtools.raw.vcf -o {sample}.concordance.raw.vcf""".format(
           gatk4=gatk4,sample=sample
        )
    return cmd

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
            outfile = """perl {extractSClip} -i {sample1}.rmdup.bam --ref_genome {ref} -p {sample1}
             perl {crest} -f {sample1}.cover -d {sample1}.rmdup.bam --ref_genome {ref} -t {ref}.2bit -p {sample1}""".format(
                extractSClip=extractSClip, sample1=sample1,crest=crest,ref=ref
            )
    return outfile
def ngs_cnv(sample1, sample2, ref, outdir, tool='control-freec',species='human',omic_type="WGS",rm_germline="false"):
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
                    config_wgs_add_control=make_freec_config.modify(sample_data,control_data,outdir,'human','WGS','Y')
                    outfile = """{control_freec} -conf {config_wgs_add_control}
                    """.format(control_freec=control_freec,config_wgs_add_control=config_wgs_add_control)
                elif rm_germline == "false":
                    sample_data = sample1+'rmdup.bam'
                    config_wgs_no_control=make_freec_config.modify(sample_data,'',outdir,'human','WGS','N')
                    outfile = """{control_freec} -conf {config_wgs_no_control}
                    """.format(control_freec=control_freec,config_wgs_no_control=config_wgs_no_control)
            elif omic_type == "WES":
                if rm_germline == "true":
                    sample_data = sample1 + 'rmdup.bam'
                    control_data = sample2 + 'rmdup.bam'
                    config_wes_add_control = make_freec_config.modify(sample_data,control_data,outdir,'human','WES','Y')
                    outfile = """{control_freec} -conf {config_wes_add_control}
                    """.format(control_freec=control_freec,config_wes_add_control=config_wes_add_control)
                elif rm_germline == "false":
                    sample_data = sample1 + 'rmdup.bam'
                    config_wes_no_control = make_freec_config.modify(sample_data, '', outdir,'human', 'WES', 'N')
                    outfile = """{control_freec} -conf {config_wes_no_control}
                    """.format(control_freec=control_freec,config_wes_no_control=config_wes_no_control)
        elif species == 'non-human':
            pass

    elif tool == 'ExomeCNV': ## WES
        pass