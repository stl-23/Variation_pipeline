### Do VQSR or Hard-filtering for GATK4 pipelines
import os
from vartools import getmyconfig
gatk4 = getmyconfig.getConfig('Variation', 'gatk4')

def vqsr(ref,vcf,vqsr_dir,out):
    vqsr_config = os.listdir(os.path.abspath(vqsr_dir) + '/vqsr_config.txt')  ### resource data for training
    with open(vqsr_config) as fh:
        for lines in fh:
            if lines.startswith('#'):
                continue
            if lines.startswith('SNP'):
                snp,args,snp_file = lines.strip().split(' ')
                new_snp_file = os.path.join(vqsr_dir,snp_file)
                snp_resource = args+' '+new_snp_file+' \\'
            elif lines.startswith('INDEL'):
                indel,args,indel_file = lines.strip.split(' ')
                new_indel_file = os.path.join(vqsr_dir,indel_file)
                indel_resource = args+' '+new_indel_file+' \\'

    cmd = """{gatk4} VariantRecalibrator -R {ref} -V {vcf} \\
{snp_resource} 
-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \\
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \\
-mode SNP \\
-rscriptFile {out}.snp.plots.R \\
--tranches-file {out}.snp.tranches \\
-O {out}.snp.recal
{gatk4} ApplyVQSR -R {ref} -V {vcf} \\
--ts_filter_level 99.0 \\
--tranches-file {out}.snp.tranches \\
-recalFile {out}.snp.recal \\
-mode SNP \\
-O {out}.snps.VQSR.gatk.vcf.gz

{gatk4} VariantRecalibrator -R {ref} -V {out}.snps.VQSR.vcf.gz \\
{indel_resource} 
-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \\
-mode INDEL --max-gaussians 6 \\
-rscriptFile {out}.snps.indel.plots.R \\
--tranches-file {out}.snps.indel.tranches \\
-O {out}.snps.indel.recal
{gatk4} ApplyVQSR -R {ref} -V {out}.snps.VQSR.vcf.gz \\
--ts_filter_level 99.0 \\
--tranches-file {out}.snps.indel.tranches \\
-recalFile {out}.snps.indel.recal \\
-mode INDEL \\
-O {out}.indel.VQSR.gatk.vcf.gz
""".format(gatk4=gatk4,ref=ref,vcf=vcf,snp_resource=snp_resource,indel_resource=indel_resource,out=out)
    return cmd

def hard_filter(ref,vcf,out):
    cmd = """{gatk4} SelectVariants -R {ref} -V {vcf} --select-type-to-include SNP -O {out}.raw.snp.vcf
{gatk4} SelectVariants -R {ref} -V {vcf} --select-type-to-include INDEL -O {out}.raw.indel.vcf
{gatk4} VariantFiltration -R {ref} -V {out}.raw.snp.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "SNP_FILTER" -O {out}.snp.genotype.gatk.vcf
{gatk4} VariantFiltration -R {ref} -V {out}.raw.indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -20.0" --filter-name "INDEL_FILTER" -O {out}.indel.genotype.gatk.vcf
{gatk4} MergeVcfs -I {out}.snp.genotype.gatk.vcf -I {out}.indel.genotype.gatk.vcf -O {out}.snp.indel.genotype.gatk.vcf
{gatk4} SelectVariants -R {ref} -V {out}.snp.indel.genotype.gatk.vcf -select "vc.isNotFiltered()" -O {out}.snp.indel.genotype.selected.gatk.vcf
""".format(gatk4=gatk4,ref=ref,vcf=vcf,out=out)
    return cmd
