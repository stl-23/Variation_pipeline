import vartools.getmyconfig as getmyconfig
gatk4 = getmyconfig.getConfig('Variation', 'gatk4')
bcftools = getmyconfig.getConfig('Variation', 'bcftools')

def mutect2(input_dir,output_dir,tmp_dir,ref,tumor,normal,interval,pon=None,germline=None,af=0.001,mem=3,threads=3,
            normals_pon=None):
    cmd_somatic = ''
    cmd_create_pon = []
    if not pon and normals_pon:
        normals_pon = [sample for sample in normals_pon.strip().split(',')]
        pon_vcfs = []
        for sample in normals_pon:
            cmd_create_pon.append("""##  Run Mutect2 in tumor-only mode for each normal sample
{gatk4} --java-options "-Xmx{mem}g" Mutect2 \\
-R {ref} \\
-I {input}/{sample}.rmdup.bam \\
-max-mnp-distance 0 \\
-O {tmp_dir}/{sample}.make.pon.vcf.gz""".format(gatk4=gatk4,mem=mem,ref=ref,sample=sample,interval=interval,
            input=input_dir,tmp_dir=tmp_dir))
            pon_vcfs.append(tmp_dir+'/'+sample+'.make.pon.vcf.gz')
            pon_vcfs_new = ' '.join(['-V '+ i for i in pon_vcfs])
        if germline:
            cmd_somatic = """### create PON ####
{gatk4} --java-options "-Xmx{mem}g" GenomicsDBImport \\
-R {ref} -L {interval} \\
{pon_vcfs_new} \\
--reader-threads {threads} \\
--genomicsdb-workspace-path {tmp_dir}/pon_db

{gatk4} --java-options "-Xmx{mem}g" CreateSomaticPanelOfNormals \\
-R {ref}  \\
--germline-resource {germline} \\
-V gendb://{tmp_dir}/pon_db \\
-O {tmp_dir}/all_sample.pon.vcf.gz 

###somatic muation###
{gatk4} --java-options "-Xmx{mem}g" Mutect2 \\
-R {ref} -L {interval} \\
-I {input}/{tumor}.rmdup.bam \\
-I {input}/{normal}.rmdup.bam \\
-tumor {tumor} \\
-normal {normal} \\
-pon {tmp_dir}/all_sample.pon.vcf.gz \\
--germline-resource {germline} \\
--af-of-alleles-not-in-resource {af} \\
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
-O {out}/{tumor}_{normal}_somatic_m2.vcf.gz

### check cross-sample contamination ###
{gatk4} --java-options "-Xmx{mem}g" SelectVariants \\
-R {ref} \\
-V {germline} \\
--select-type-to-include SNP \\
--restrict-alleles-to BIALLELIC \\
-O {out}/germline.SNP_biallelic.vcf.gz
{gatk4} --java-options "-Xmx{mem}g" GetPileupSummaries \\
-I {input}/{tumor}.rmdup.bam \\
-L {interval} \\
-V {out}/germline.SNP_biallelic.vcf.gz \\
-O {out}/{tumor}.pileups.table
{gatk4} --java-options "-Xmx{mem}g" GetPileupSummaries \\
-I {input}/{normal}.rmdup.bam \\
-L {interval} \\
-V {out}/germline.SNP_biallelic.vcf.gz
-O {out}/{normal}.pileups.table
{gatk4} --java-options "-Xmx{mem}g" CalculateContamination \\
-I {out}/{tumor}.pileups.table \\
-matched {out}/{normal}.pileups.table \\
-O {out}/{tumor}.calculatecontamination.table

### filter
{gatk4} --java-options "-Xmx{mem}g" FilterMutectCalls \\
-R {ref} \\
-V {out}/{tumor}_{normal}_somatic_m2.vcf.gz \\
--contamination-table {out}/{tumor}.calculatecontamination.table \\
-O {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz

zcat {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz |egrep '#|PASS' | {bcftools} view -Oz -o {out}/{tumor}_{normal}_somatic_m2.filter.pass.vcf.gz
{gatk4} IndexFeatureFile -I  {out}/{tumor}_{normal}_somatic_m2.filter.pass.vcf.gz
""".format(gatk4=gatk4,mem=mem,ref=ref,interval=interval,pon_vcfs_new=pon_vcfs_new,out=output_dir,
           input=input_dir,tumor=tumor,normal=normal,germline=germline,af=af,bcftools=bcftools,
           tmp_dir=tmp_dir,threads=threads)
        else:
            cmd_somatic = """### create PON ###
{gatk4} --java-options "-Xmx{mem}g" GenomicsDBImport \\
-R {ref} -L {interval} \\
{pon_vcfs_new} \\
--reader-threads {threads} \\
--genomicsdb-workspace-path {tmp_dir}/pon_db

{gatk4} --java-options "-Xmx{mem}g" CreateSomaticPanelOfNormals \\
-R {ref}  \\
-V gendb://{tmp_dir}/pon_db \\
-O {tmp_dir}/all_sample.pon.vcf.gz 

### somatic mutation ###
{gatk4} --java-options "-Xmx{mem}g"  Mutect2 \\
-R {ref} -L {interval} \\
-I {input}/{tumor}.rmdup.bam \\
-I {input}/{normal}.rmdup.bam \\
-tumor {tumor} \\
-normal {normal} \\
-pon {tmp_dir}/all_sample.pon.vcf.gz \\
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
-O {out}/{tumor}_{normal}_somatic_m2.vcf.gz

### filter ###
{gatk4} --java-options "-Xmx{mem}g" FilterMutectCalls \\
-R {ref} \\
-V {out}/{tumor}_{normal}_somatic_m2.vcf.gz \\
-O {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz
zcat {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz |egrep '#|PASS' | {bcftools} view -Oz -o {out}/{tumor}_{normal}_somatic_m2.filter.pass.vcf.gz
{gatk4} IndexFeatureFile -I  {out}/{tumor}_{normal}_somatic_m2.filter.pass.vcf.gz
""".format(gatk4=gatk4, mem=mem, ref=ref, interval=interval, pon_vcfs_new=pon_vcfs_new, out=output_dir,
           input=input_dir, tumor=tumor, normal=normal,bcftools=bcftools,tmp_dir=tmp_dir,threads=threads)
    elif pon and not normals_pon:
        if germline:
            cmd_somatic = """### somatic mutation ###
{gatk4} --java-options "-Xmx{mem}g"  Mutect2 \\
-R {ref} -L {interval} \\
-I {input}/{tumor}.rmdup.bam \\
-I {input}/{normal}.rmdup.bam \\
-tumor {tumor} \\
-normal {normal} \\
-pon {pon} \\
--germline-resource {germline} \\
--af-of-alleles-not-in-resource {af} \\
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
-O {out}/{tumor}_{normal}_somatic_m2.vcf.gz

### check cross-sample contamination ### 
{gatk4} --java-options "-Xmx{mem}g" SelectVariants \\
-R {ref} \\
-V {germline} \\
--select-type-to-include SNP \\
--restrict-alleles-to BIALLELIC \\
-O {out}/germline.SNP_biallelic.vcf.gz
{gatk4} --java-options "-Xmx{mem}g" GetPileupSummaries \\
-I {input}/{tumor}.rmdup.bam \\
-L {interval} \\
-V {out}/germline.SNP_biallelic.vcf.gz \\
-O {out}/{tumor}.pileups.table
{gatk4} --java-options "-Xmx{mem}g" GetPileupSummaries \\
-I {input}/{normal}.rmdup.bam \\
-L {interval} \\
-V {out}/germline.SNP_biallelic.vcf.gz \\
-O {out}/{normal}.pileups.table
{gatk4} --java-options "-Xmx{mem}g" CalculateContamination \\
-I {out}/{tumor}.pileups.table \\
-matched {out}/{normal}.pileups.table \\
-O {out}/{tumor}.calculatecontamination.table

### filter
{gatk4} --java-options "-Xmx{mem}g" FilterMutectCalls \\
-R {ref} \\
-V {out}/{tumor}_{normal}_somatic_m2.vcf.gz \\
--contamination-table {out}/{tumor}.calculatecontamination.table \\
-O {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz
zcat {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz |egrep '#|PASS' | {bcftools} view -Oz -o {out}/{tumor}_{normal}_somatic_m2.filter.pass.vcf.gz
{gatk4} IndexFeatureFile -I  {out}/{tumor}_{normal}_somatic_m2.filter.pass.vcf.gz
""".format(gatk4=gatk4,mem=mem,ref=ref,interval=interval,out=output_dir,pon=pon,
           input=input_dir,tumor=tumor,normal=normal,germline=germline,af=af,bcftools=bcftools)
        else:
            cmd_somatic = """### somatic mutation ###
{gatk4} --java-options "-Xmx{mem}g"  Mutect2 \\
-R {ref} -L {interval} \\
-I {input}/{tumor}.rmdup.bam \\
-I {input}/{normal}.rmdup.bam \\
-tumor {tumor} \\
-normal {normal} \\
-pon {pon} \\
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
-O {out}/{tumor}_{normal}_somatic_m2.vcf.gz
###filter###
{gatk4} --java-options "-Xmx{mem}g" FilterMutectCalls \\
-R {ref} \\
-V {out}/{tumor}_{normal}_somatic_m2.vcf.gz \\
-O {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz
zcat {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz |egrep '#|PASS' | {bcftools} view -Oz -o {out}/{tumor}_{normal}_somatic_m2.filter.pass.vcf.gz
{gatk4} IndexFeatureFile -I  {out}/{tumor}_{normal}_somatic_m2.filter.pass.vcf.gz
""".format(gatk4=gatk4, mem=mem, ref=ref, interval=interval, out=output_dir, pon=pon,
                       input=input_dir, tumor=tumor, normal=normal,bcftools=bcftools)
    elif not pon and not normals_pon:
        if germline:
            cmd_somatic = """### somatic mutation ###
{gatk4} --java-options "-Xmx{mem}g"  Mutect2 \\
-R {ref} -L {interval} \\
-I {input}/{tumor}.rmdup.bam \\
-I {input}/{normal}.rmdup.bam \\
-tumor {tumor} \\
-normal {normal} \\
--germline-resource {germline} \\
--af-of-alleles-not-in-resource {af} \\
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
-O {out}/{tumor}_{normal}_somatic_m2.vcf.gz

### check cross-sample contamination ### 
{gatk4} --java-options "-Xmx{mem}g" SelectVariants \\
-R {ref} \\
-V {germline} \\
--select-type-to-include SNP \\
--restrict-alleles-to BIALLELIC \\
-O {out}/germline.SNP_biallelic.vcf.gz
{gatk4} --java-options "-Xmx{mem}g" GetPileupSummaries \\
-I {input}/{tumor}.rmdup.bam \\
-L {interval} \\
-V {out}/germline.SNP_biallelic.vcf.gz
-O {out}/{tumor}.pileups.table
{gatk4} --java-options "-Xmx{mem}g" GetPileupSummaries \\
-I {input}/{normal}.rmdup.bam \\
-L {interval} \\
-V {out}/germline.SNP_biallelic.vcf.gz \\
-O {out}/{normal}.pileups.table
{gatk4} --java-options "-Xmx{mem}g" CalculateContamination \\
-I {out}/{tumor}.pileups.table \\
-matched {out}/{normal}.pileups.table \\
-O {out}/{tumor}.calculatecontamination.table

### filter
{gatk4} --java-options "-Xmx{mem}g" FilterMutectCalls \\
-R {ref} \\
-V {out}/{tumor}_{normal}_somatic_m2.vcf.gz \\
--contamination-table {out}/{tumor}.calculatecontamination.table \\
-O {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz
zcat {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz |egrep '#|PASS' | {bcftools} view -Oz -o {out}/{tumor}_{normal}_somatic_m2.filter.pass.vcf.gz
{gatk4} IndexFeatureFile -I  {out}/{tumor}_{normal}_somatic_m2.filter.pass.vcf.gz
""".format(gatk4=gatk4,mem=mem,ref=ref,interval=interval,out=output_dir,
           input=input_dir,tumor=tumor,normal=normal,germline=germline,af=af,bcftools=bcftools)
        else:
            cmd_somatic = """### somatic mutation ###
{gatk4} --java-options "-Xmx{mem}g"  Mutect2 \\
-R {ref} -L {interval} \\
-I {input}/{tumor}.rmdup.bam \\
-I {input}/{normal}.rmdup.bam \\
-tumor {tumor} \\
-normal {normal} \\
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
-O {out}/{tumor}_{normal}_somatic_m2.vcf.gz
###filter###
{gatk4} --java-options "-Xmx{mem}g" FilterMutectCalls \\
-R {ref} \\
-V {out}/{tumor}_{normal}_somatic_m2.vcf.gz \\
-O {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz
zcat {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz |egrep '#|PASS' | {bcftools} view -Oz -o {out}/{tumor}_{normal}_somatic_m2.filter.pass.vcf.gz
{gatk4} IndexFeatureFile -I  {out}/{tumor}_{normal}_somatic_m2.filter.pass.vcf.gz
""".format(gatk4=gatk4, mem=mem, ref=ref, interval=interval, out=output_dir,
                       input=input_dir, tumor=tumor, normal=normal,bcftools=bcftools)

    return cmd_create_pon,cmd_somatic

#def varscan():


#def strelka():