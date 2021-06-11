import os
from vartools import getmyconfig
annovar_dir = getmyconfig.getConfig('Variation', 'ANNOVAR')
gff3ToGenePred = getmyconfig.getConfig('Variation','gff3ToGenePred')
hg19_db = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), 'database/genomicsdb/hg19/annovar')
hg38_db = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), 'database/genomicsdb/hg38/annovar')

def annotation(tool, ref, vcf, gff3, out, species='hg19'):
    cmd = ''
    if tool == 'annovar':
        if species == 'hg19':
            cmd = """perl {annovar_dir}/convert2annovar.pl -format vcf4 {vcf} > {out}.avinput
perl {annovar_dir}/annotate_variation.pl -buildver hg19 -geneanno -dbtype refGene {out}.avinput {humandb} --outfile {out}
""".format(annovar_dir=annovar_dir,vcf=vcf,out=out,humandb=hg19_db)
        elif species == 'hg38':
            cmd = """perl {annovar_dir}/convert2annovar.pl -format vcf4 {vcf} > {out}.avinput
perl {annovar_dir}/annotate_variation.pl -buildver hg38 -geneanno -dbtype refGene {out}.avinput {humandb} --outfile {out}
""".format(annovar_dir=annovar_dir,vcf=vcf,out=out,humandb=hg38_db)
        else:
            cmd = """perl {annovar_dir}/convert2annovar.pl -format vcf4 {vcf} > {out}.avinput
{gff3ToGenePred} {gff3} {species}/{species}_refGene.txt
perl {annovar_dir}/retrieve_seq_from_fasta.pl --format refGene --seqfile {ref} {species}/{species}_refGene.txt --out {species}/{species}_efGeneMrna.fa
perl {annovar_dir}/annotate_variation.pl -dbtype refGene {out}.avinput {species} --outfile {out}
""".format(annovar_dir=annovar_dir,vcf=vcf,ref=ref,out=out,gff3ToGenePred=gff3ToGenePred,gff3=gff3,species=species)
    return cmd