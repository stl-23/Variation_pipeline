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
splitSNPindelVCF = getmyconfig.getConfig('Variation','splitSNPindelVCF')
makeGraph = getmyconfig.getConfig('Variation','makeGraph')

## known_site='--known-sites /path/to/ref1.vcf --known-sites /path/to/ref2.vcf ....'
def snp_indel_samtools(ref, input, sample, v_valling, bcftools_filter):
    outfile = ''
    #samtools_p = 'mpileup -C 50 -m 2 -F 0.002 -d 1000 -u -f'
    # vcfutils_p = 'varFilter -Q 20 -d 4 -D 1000'
    bcftools_mpileup = 'mpileup -d 1000 -Ov -f'
    bcftools_call = 'call -mv -Oz -o'
    ### Hard filtering
    if v_valling == 'single':  ## input is a single bam file
        outfile = """{bcftools} {bcftools_mpileup} {ref} {input_bam} | {bcftools} {bcftools_call} {sample}.all.vcf.gz
{bcftools} filter {bcftools_filter} {sample}.all.vcf.gz -o {sample}.filter.vcf.gz
python {splitSNPindelVCF} {sample}.filter.vcf.gz {sample}.samtools
{bcftools} view {sample}.samtools.snp.vcf -Oz -o {sample}.samtools.snp.vcf.gz
{bcftools} view {sample}.samtools.indel.vcf -Oz -o {sample}.samtools.indel.vcf.gz
{gatk4} IndexFeatureFile -I {sample}.samtools.snp.vcf.gz 
{gatk4} IndexFeatureFile -I {sample}.samtools.indel.vcf.gz
rm {sample}.samtools.snp.vcf {sample}.samtools.indel.vcf 
        """.format(bcftools=bcftools,bcftools_mpileup=bcftools_mpileup,bcftools_call=bcftools_call,
                   ref=ref,sample=sample,input_bam=input, gatk4=gatk4,
                   bcftools_filter=bcftools_filter,splitSNPindelVCF=splitSNPindelVCF)
    elif v_valling == 'join': ## input is a list of bam files
        input_bams = ' '.join(input)
        outfile = """{bcftools} {bcftools_mpileup} {ref} {input_bams} | {bcftools} {bcftools_call} {sample}.all.vcf.gz
{bcftools} filter {bcftools_filter} {sample}.all.vcf.gz -o {sample}.filter.vcf.gz
python {splitSNPindelVCF} {sample}.filter.vcf.gz {sample}.samtools
{bcftools} view {sample}.samtools.snp.vcf -Oz -o {sample}.samtools.snp.vcf.gz
{bcftools} view {sample}.samtools.indel.vcf -Oz -o {sample}.samtools.indel.vcf.gz
{gatk4} IndexFeatureFile -I {sample}.samtools.snp.vcf.gz 
{gatk4} IndexFeatureFile -I {sample}.samtools.indel.vcf.gz
rm {sample}.samtools.snp.vcf {sample}.samtools.indel.vcf """.format(
            bcftools=bcftools, bcftools_mpileup=bcftools_mpileup, bcftools_call=bcftools_call,
            ref=ref, sample=sample, input_bams=input_bams,gatk4=gatk4,
            bcftools_filter=bcftools_filter, splitSNPindelVCF=splitSNPindelVCF
        )
    return outfile

def snp_indel_gatk(ref, input, sample, gvcf, bqsr_dir):
    outfile = ''
    if bqsr_dir:  #### BQSR (Recalibration Base Quality Score)
        lst = os.listdir(os.path.abspath(bqsr_dir))
        known_site = ' '.join(["--known-sites "+bqsr_dir+"/"+vcf for vcf in lst if vcf.endswith('.vcf.gz')])
        if gvcf == "T":
            outfile = """{gatk4} BaseRecalibrator -R {ref} -I {input} {known_site} -O {sample}.recal.table
{gatk4} ApplyBQSR --bqsr-recal-file {sample}.recal.table -R {ref} -I {input} -O {sample}.BQSR.bam
{gatk4} HaplotypeCaller -R {ref} -I {sample}.BQSR.bam -ERC GVCF -O {sample}.g.vcf
                    """.format(gatk4=gatk4, ref=ref, input=input, sample=sample, known_site=known_site)

        elif gvcf == "F":
            outfile = """{gatk4} BaseRecalibrator -R {ref} -I {input} {known_site} -O {sample}.recal.table
{gatk4} ApplyBQSR --bqsr-recal-file {sample}.recal.table -R {ref} -I {input} -O {sample}.BQSR.bam
{gatk4} HaplotypeCaller -R {ref} -I {sample}.BQSR.bam -O {sample}.vcf
                        """.format(gatk4=gatk4, ref=ref, input=input, sample=sample, known_site=known_site)
    elif not bqsr_dir:
        if gvcf == "T":
            outfile = """{gatk4} HaplotypeCaller -R {ref} -I {input} -ERC GVCF -O {sample}.g.vcf""".format(gatk4=gatk4, ref=ref, input=input, sample=sample)
        elif gvcf == "F":
            outfile = """{gatk4} HaplotypeCaller -R {ref} -I {input} -O {sample}.vcf""".format(gatk4=gatk4, ref=ref, input=input, sample=sample)
    return outfile

def samtool_gatk_combine(sample):  ## Samtools and GATK pipelines
    cmd = """{gatk4} SelectVariants --variant {sample}.snps.gatk.vcf.gz --concordance {sample}.samtools.snp.vcf.gz -O {sample}.concordance.snp.vcf.gz
{gatk4} SelectVariants --variant {sample}.indel.gatk.vcf.gz --concordance {sample}.samtools.indel.vcf.gz -O {sample}.concordance.indel.vcf.gz""".format(
        gatk4=gatk4,sample=sample)
    return cmd

def ngs_sv(sample1,sample2,ref,out,tool='breakdancer',rm_germline="false"): ## sample1:disease/somatic; sample2:control
    outfile = ''
    breakdancer_p = '-q 20 -d'
    input_path = os.path.dirname(sample1)
    sample1_name = os.path.basename(sample1).replace('.rmdup.bam', '')
    #sample2_name = os.path.basename(sample2).replace('.rmdup.bam', '')
    #crest_p = ''
    if tool == 'breakdancer':  ## Only use paired end reads
        if rm_germline == "true": ## somatic SV
            pass
        elif rm_germline == "false": ## germline/common SV
            outfile = """perl {bam2cfg} -q 20 -c 4 -g {sample1} > {sample1}.cfg
{breakdancer} {breakdancer_p} {sample1_name} {sample1}.cfg > {out}.raw.ctx""".format(
                bam2cfg=bam2cfg,sample1=sample1,breakdancer=breakdancer,breakdancer_p=breakdancer_p,
                out=out,sample1_name=sample1_name
            )
    elif tool == 'crest':      ## Can use both paired end and singe end reads
        if rm_germline == "true":  ## somatic SV
            outfile = """ln -s {input}/{sample1_name}.rmdup.bai {input}/{sample1_name}.rmdup.bam.bai
perl {extractSClip} -i {sample1} --ref_genome {ref} -o {out}
perl {crest} -f {sample1}.cover -d {sample1} -g {sample2} --ref_genome {ref} -t {ref}.2bit -o {out}
            """.format(
                extractSClip=extractSClip,sample1=sample1,sample2=sample2,
                crest=crest,ref=ref,sample1_name=sample1_name,out=out,input=input_path
            )
        elif rm_germline == "false": ## germline/common SV
            outfile = """ln -s {input}/{sample1_name}.rmdup.bai {input}/{sample1_name}.rmdup.bam.bai
perl {extractSClip} -i {sample1} --ref_genome {ref} -o {out}
perl {crest} -f {sample1}.cover -d {sample1} --ref_genome {ref} -t {ref}.2bit -o {out}""".format(
                extractSClip=extractSClip, sample1=sample1,crest=crest,ref=ref,sample1_name=sample1_name,out=out,
                input=input_path
            )
    return outfile
def ngs_cnv(sample1, sample2, ref, outdir, tool='control-freec',species='human',stragety="WGS",rm_germline="false"):
    outfile = ''
    sample1_name = os.path.basename(sample1).replace('.rmdup.bam', '')
    out_name = outdir+'/'+sample1_name
    if tool == 'cnvnator':
        bin_size = '100'  ## sequencing depth=20-30x: bin size = 100; 2-3x: 500; 100x:30
        ref_dir = os.path.dirname(ref)
        outfile = """{cnvnator} -root {out_name}.root -tree {sample1} 
## {cnvnator} -root {out_name}.root -eval {bin_size}  # if bin size was suitable, the value of mean to segma would be 4-5
{cnvnator} -root {out_name}.root -his {bin_size} -d {ref_dir}/split_fa
{cnvnator} -root {out_name}.root -stat {bin_size}
{cnvnator} -root {out_name}.root -partition {bin_size}
{cnvnator} -root {out_name}.root -call {bin_size} > {out_name}.cnv.call.txt
perl {cnvnator2VCF} {out_name}.cnv.call.txt > {out_name}.cnv.vcf
        """.format(cnvnator=cnvnator,out_name=out_name,sample1=sample1,ref_dir=ref_dir,cnvnator2VCF=cnvnator2VCF,
                   bin_size=bin_size)
        return outfile
    elif tool == 'control-freec':
        if species == 'human':
            if stragety == "WGS":
                if rm_germline == "true":
                    sample_data = sample1+'.rmdup.bam'
                    control_data = sample2+'.rmdup.bam'
                    config_wgs_add_control=make_freec_config.modify(sample_data,control_data,ref,outdir,'human','WGS','Y')
                    outfile = """{control_freec} -conf ./config.list
awk '$3!=-1 ' {sample_data}_ratio.txt > {sample_data}_ratio_noNA.txt
cat {makeGraph} | R --slave --args 2 {sample_data}_ratio_noNA.txt
                    """.format(control_freec=control_freec,sample_data=sample_data,makeGraph=makeGraph)
                    return config_wgs_add_control, outfile
                elif rm_germline == "false":
                    pre = sample1.rstrip('/').split('/')[-1]
                    config_wgs_no_control=make_freec_config.modify(sample1,'',ref,outdir,'human','WGS','N')
                    outfile = """{control_freec} -conf {pre}_config_wgs_no_control.list
awk '$3!=-1 ' {sample_data}_ratio.txt > {sample_data}_ratio_noNA.txt
cat {makeGraph} | R --slave --args 2 {sample_data}_ratio_noNA.txt
                    """.format(control_freec=control_freec,pre=pre,sample_data=sample_data,makeGraph=makeGraph)
                    return config_wgs_no_control, outfile
            elif stragety == "WES":
                if rm_germline == "true":
                    sample_data = sample1 + '.rmdup.bam'
                    control_data = sample2 + '.rmdup.bam'
                    config_wes_add_control = make_freec_config.modify(sample_data,control_data,ref,outdir,'human','WES','Y')
                    outfile = """{control_freec} -conf ./config.list
awk '$3!=-1 ' {sample_data}_ratio.txt > {sample_data}_ratio_noNA.txt
cat {makeGraph} | R --slave --args 2 {sample_data}_ratio_noNA.txt
                    """.format(control_freec=control_freec,sample_data=sample_data,makeGraph=makeGraph)
                    return config_wes_add_control, outfile
                elif rm_germline == "false":
                    pre = sample1.rstrip('/').split('/')[-1]
                    config_wes_no_control = make_freec_config.modify(sample1, '', ref,outdir,'human', 'WES', 'N')
                    outfile = """{control_freec} -conf {pre}_config_wgs_no_control.list
awk '$3!=-1 ' {sample_data}_ratio.txt > {sample_data}_ratio_noNA.txt
cat {makeGraph} | R --slave --args 2 {sample_data}_ratio_noNA.txt
                    """.format(control_freec=control_freec,pre=pre,sample_data=sample_data,makeGraph=makeGraph)
                    return config_wes_no_control, outfile
        elif species == 'non-human':
            pass

    elif tool == 'ExomeCNV': ## WES
        pass