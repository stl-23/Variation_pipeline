import vartools.getmyconfig as getmyconfig

bcftools = getmyconfig.getConfig('Variation', 'bcftools')
gatk4 = getmyconfig.getConfig('Variation', 'gatk4')

def merge(files,type,prefix,*gvcf_p):
    cmd = ''
    if type == 'vcf':
        combine_vcf = ' '.join([i + 'sort.gz' for i in files])
        for vcf in files:
            cmd = """{bcftools} view {vcf} -Oz -o {vcf}.gz
{bcftools} sort {vcf}.gz -o {vcf}.sort.gz
{bcftools} index {vcf}.sort.gz """.format(
                bcftools=bcftools,vcf=vcf,
            )
        cmd+="""{bcftools} merge {combine_vcf} -o {prefix}.merged.vcf""".format(
            bcftools=bcftools, combine_vcf=combine_vcf,prefix=prefix
        )
    elif type == 'gvcf':
        combine_gvcf = ' '.join(['-V '+i for i in files])
        mem, ref, genomicsdb, chr_list, b_size, map_file, reader_threads, num, tmp = gvcf_p[:]

        if len(files) < 1000:
            cmd = """{gatk4} --java-options "-Xmx{mem}g" CombineGVCFs -R {ref} {combine_gvcf} -O {prefix}.combined.g.vcf
{gatk4} --java-options "-Xmx{mem}g" GenotypeGVCFs -R {ref} -V {prefix}.combined.g.vcf -G StandardAnnotation -new-qual \\
-O {prefix}.combined.vcf""".format(gatk4=gatk4,mem=mem,ref=ref,combine_gvcf=combine_gvcf, prefix=prefix)

        elif len(files) >= 1000:
            cmd = """{gatk4} --java-options "-Xmx{mem}g" GenomicsDBImport --genomicsdb-workspace-path {genomicsdb} \\
--intervals {chr_list} --batch-size {b_size} --sample-name-map {map_file} --reader-threads {reader_threads} \\
--max-num-intervals-to-import-in-parallel {num} --tmp-dir {tmp}
{gatk4} --java-options "-Xmx{mem}g" GenotypeGVCFs -R {ref} -V gendb:{genomicsdb} -G StandardAnnotation -new-qual \\
-O {prefix}.combined.vcf""".format(gatk4=gatk4,mem=mem,genomicsdb=genomicsdb,chr_list=chr_list,b_size=b_size,
                                   map_file=map_file,reader_threads=reader_threads,num=num,tmp=tmp,
                                   ref=ref,prefix=prefix)
    return cmd
