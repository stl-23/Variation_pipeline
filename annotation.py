import os
from tools import getmyconfig
annovar_dir = getmyconfig.getConfig('Variation', 'ANNOVAR')
gff3ToGenePred = getmyconfig.getConfig('Variation','gff3ToGenePred')
hg19_db = '/tool/annovar/humandb/'

def annotation(tool='annovar', vcf, ref, gff3, species='hg19', out):
    cmd = ''
    if tool == 'annovar':
        if species == 'hg19':
            cmd = """perl {annovar_dir}/convert2annovar.pl -format vcf4 {vcf} > {sample}.avinput
perl {annovar_dir}/annotate_variation.pl -buildver hg19 -geneanno -dbtype refGene {sample}.avinput {humandb} --outfile {sample}
""".format(annovar_dir=annovar_dir,vcf=vcf,sample=out,humandb=hg19_db)

        else:
            cmd = """perl {annovar_dir}/convert2annovar.pl -format vcf4 {vcf} > {sample}.avinput
{gff3ToGenePred} {gff3} {species}/{species}_refGene.txt
perl {annovar_dir}/retrieve_seq_from_fasta.pl --format refGene --seqfile {ref} {species}/{species}_refGene.txt --out {species}/{species}_efGeneMrna.fa
perl {annovar_dir}/annotate_variation.pl -dbtype refGene {sample}.avinput {species} --outfile {sample}
""".format(annovar_dir=annovar_dir,vcf=vcf,ref=ref,sample=out,gff3ToGenePred=gff3ToGenePred,gff3=gff3,species=species)
    return cmd