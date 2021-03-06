### Do VQSR or Hard-filtering for GATK4 pipelines
import os
import vartools.getmyconfig as getmyconfig
gatk4 = getmyconfig.getConfig('Variation', 'gatk4')

def vqsr(ref,vcf,vqsr_dir,out):
    vqsr_config = os.path.join(os.path.abspath(vqsr_dir), 'vqsr_config.txt')  ### resource data for training
    with open(vqsr_config) as fh:
        snp_resource = []
        indel_resource = []
        for lines in fh:
            if lines.startswith('#'):
                continue
            if lines.startswith('SNP'):
                snp,args,snp_file = lines.strip().split(' ')
                new_snp_file = os.path.join(vqsr_dir,snp_file)
                snp_resource.append(args+' '+new_snp_file+' \\')
            elif lines.startswith('INDEL'):
                indel,args,indel_file = lines.strip().split(' ')
                new_indel_file = os.path.join(vqsr_dir,indel_file)
                indel_resource.append(args+' '+new_indel_file+' \\')
        snp_resource_all = '\n'.join(snp_resource).strip('\\')
        indel_resource_all = '\n'.join(indel_resource).strip('\\')
        cmd = """{gatk4} VariantRecalibrator -R {ref} -V {vcf} \\
{snp_resource_all} \\
--trust-all-polymorphic \\
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \\
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \\
-mode SNP --max-gaussians 6 \\
#--rscript-file {out}.snp.plots.R \\
--tranches-file {out}.snp.tranches \\
-O {out}.snp.recal
{gatk4} ApplyVQSR -R {ref} -V {vcf} \\
--truth-sensitivity-filter-level 99.7 \\
--tranches-file {out}.snp.tranches \\
--recal-file {out}.snp.recal \\
-mode SNP \\
-O {out}.snps.gatk.vcf.gz

{gatk4} VariantRecalibrator -R {ref} -V {out}.snps.gatk.vcf.gz \\
{indel_resource_all} \\
--trust-all-polymorphic \\
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \\
-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \\ 
-mode INDEL --max-gaussians 4 \\
#--rscript-file {out}.indel.plots.R \\
--tranches-file {out}.indel.tranches \\
-O {out}.indel.recal
{gatk4} ApplyVQSR -R {ref} -V {out}.snps.gatk.vcf.gz \\
--truth-sensitivity-filter-level 99.7 \\
--tranches-file {out}.indel.tranches \\
--recal-file {out}.indel.recal \\
-mode INDEL \\
-O {out}.indel.gatk.vcf.gz
""".format(gatk4=gatk4,ref=ref,vcf=vcf,snp_resource_all=snp_resource_all,indel_resource_all=indel_resource_all,out=out)
    return cmd

def hard_filter(ref,vcf,out):
    cmd = """{gatk4} SelectVariants --java-options "-Xmx8G" -R {ref} -V {vcf} --select-type-to-include SNP -O {out}.raw.snp.vcf
{gatk4} SelectVariants --java-options "-Xmx8G" -R {ref} -V {vcf} --select-type-to-include INDEL -O {out}.raw.indel.vcf
{gatk4} VariantFiltration -R {ref} -V {out}.raw.snp.vcf \\
-filter "QD < 2.0" --filter-name "QD2" \\
-filter "QUAL < 30.0" --filter-name "QUAL30" \\
-filter "FS > 60.0" --filter-name "FS60" \\
-filter "MQ < 40.0" --filter-name "MQ40" \\
-filter "SOR > 3.0" --filter-name "SOR3" \\
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
-O {out}.snp.gatk.vcf.gz
{gatk4} VariantFiltration -R {ref} -V {out}.raw.indel.vcf \\
-filter "QD < 2.0"  --filter-name "QD2" \\
-filter "QUAL < 30.0" --filter-name "QUAL30" \\
-filter "FS > 200.0" --filter-name "FS200" \\
-filter "SOR > 10.0" --filter-name "SOR10" \\
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-8" \\
-O {out}.indel.gatk.vcf.gz
#{gatk4} MergeVcfs -I {out}.snp.genotype.gatk.vcf -I {out}.indel.genotype.gatk.vcf -O {out}.snp.indel.genotype.gatk.vcf
#{gatk4} SelectVariants -R {ref} -V {out}.snp.indel.genotype.gatk.vcf -select "vc.isNotFiltered()" -O {out}.snp.indel.genotype.selected.gatk.vcf
""".format(gatk4=gatk4,ref=ref,vcf=vcf,out=out)
    return cmd
