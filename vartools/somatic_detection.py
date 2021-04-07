import vartools.getmyconfig as getmyconfig
gatk4 = getmyconfig.getConfig('Variation', 'gatk4')

def mutect2(input_dir,output_dir,ref,tumor,normal,interval,pon='',germline='',af=0.001,mem=3,*normals_pon):
    cmd_somatic = ''
    cmd_create_pon = []
    if not pon and normals_pon:
        pon_vcfs = []
        for sample in normals_pon:
            cmd_create_pon.append("""{gatk4} --java-options "-Xmx{mem}g Mutect2 \\
-R {ref} -I {input}/{sample}.rmdup.bam -tumor {sample} --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
-L {interval} -O {out}/{sample}.prepon.vcf.gz""".format(gatk4=gatk4,mem=mem,ref=ref,sample=sample,interval=interval,
            input=input_dir,out=output_dir))
            pon_vcfs.append(output_dir+'/'+sample+'.prepon.vcf.gz')
            pon_vcfs_new = ' '.join(['-V '+ i for i in pon_vcfs])
        if germline:
            cmd_somatic = """{gatk4} --java-options "-Xmx{mem}g" CreateSomaticPanelOfNormals \\
-R {ref} -L {interval} \\
{pon_vcfs_new} \\
-O {out}/all_sample.prepon.vcf.gz \\

{gatk4} --java-options "-Xmx{mem}g"  Mutect2 \\
-R {ref} -L {interval} \\
-I {input}/{tumor}.rmdup.bam \\
-I {input}/{normal}.rmdup.bam \\
-tumor {tumor} \\
-normal {normal} \\
-pon {out}/all_sample.prepon.vcf.gz \\
--germline-resource {germline} \\
--af-of-alleles-not-in-resource {af} \\
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
-O {out}/{tumor}_{normal}_somatic_m2.vcf.gz

{gatk4} --java-options "-Xmx{mem}g" FilterMutectCalls \\
-V {out}/{tumor}_{normal}_somatic_m2.vcf.gz \\
-O {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz
""".format(gatk4=gatk4,mem=mem,ref=ref,interval=interval,pon_vcfs_new=pon_vcfs_new,out=output_dir,
           input=input_dir,tumor=tumor,normal=normal,germline=germline,af=af)
        else:
            cmd_somatic = """{gatk4} --java-options "-Xmx{mem}g" CreateSomaticPanelOfNormals \\
-R {ref} -L {interval} \\
{pon_vcfs_new} \\
-O {out}/all_sample.prepon.vcf.gz \\

{gatk4} --java-options "-Xmx{mem}g"  Mutect2 \\
-R {ref} -L {interval} \\
-I {input}/{tumor}.rmdup.bam \\
-I {input}/{normal}.rmdup.bam \\
-tumor {tumor} \\
-normal {normal} \\
-pon {out}/all_sample.prepon.vcf.gz \\
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
-O {out}/{tumor}_{normal}_somatic_m2.vcf.gz

{gatk4} --java-options "-Xmx{mem}g" FilterMutectCalls \\
-V {out}/{tumor}_{normal}_somatic_m2.vcf.gz \\
-O {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz
""".format(gatk4=gatk4, mem=mem, ref=ref, interval=interval, pon_vcfs_new=pon_vcfs_new, out=output_dir,
                       input=input_dir, tumor=tumor, normal=normal)
    elif pon and not normals_pon:
        if germline:
            cmd_somatic = """{gatk4} --java-options "-Xmx{mem}g"  Mutect2 \\
-R {ref} -L {interval} \\
-I {input}/{tumor}.rmdup.bam \\
-I {input}/{normal}.rmdup.bam \\
-tumor {tumor} \\
-normal {normal} \\
-pon {input}/{pon} \\
--germline-resource {germline} \\
--af-of-alleles-not-in-resource {af} \\
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
-O {out}/{tumor}_{normal}_somatic_m2.vcf.gz

{gatk4} --java-options "-Xmx{mem}g" FilterMutectCalls \\
-V {out}/{tumor}_{normal}_somatic_m2.vcf.gz \\
-O {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz
""".format(gatk4=gatk4,mem=mem,ref=ref,interval=interval,out=output_dir,pon=pon,
           input=input_dir,tumor=tumor,normal=normal,germline=germline,af=af)
        else:
            cmd_somatic = """{gatk4} --java-options "-Xmx{mem}g"  Mutect2 \\
-R {ref} -L {interval} \\
-I {input}/{tumor}.rmdup.bam \\
-I {input}/{normal}.rmdup.bam \\
-tumor {tumor} \\
-normal {normal} \\
-pon {input}/{pon} \\
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
-O {out}/{tumor}_{normal}_somatic_m2.vcf.gz

{gatk4} --java-options "-Xmx{mem}g" FilterMutectCalls \\
-V {out}/{tumor}_{normal}_somatic_m2.vcf.gz \\
-O {out}/{tumor}_{normal}_somatic_m2.filter.vcf.gz
""".format(gatk4=gatk4, mem=mem, ref=ref, interval=interval, out=output_dir, pon=pon,
                       input=input_dir, tumor=tumor, normal=normal)

    return cmd_create_pon,cmd_somatic

#def varscan():


#def strelka():