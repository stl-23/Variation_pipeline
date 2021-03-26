#!/usr/bin/env python
import os
import sys
import argparse
import mapping
import annotation
import ngs_vars
import tgs_vars
from vartools import parsering, getmyconfig

script_path = os.path.abspath("./")
make_index = os.path.join(script_path, '/tools/index.sh')

def main():
    input_path = os.path.abspath(self.inputs) + '/'
    out_path = os.path.abspath(self.outputs) + '/'
    # lst = os.listdir(input_path)
    outfile = []
    samples, input1, input2, output1, output2 = parsering.parse_short_read_dir(self.input_path, self.out_path,
                                                                               self.type)

    if platform == "NGS":
        out1,out2,out3 = ngs_vars.ngs_snp_indel(ref, sample, tool, combine,BQSR,gcvf,known_site)
        if combine == "N" and tool == 'samtools':
            fw = open('samtools_snp_indel.sh','w').write(out1)
        elif combine == 'N' and tool == 'gatk':
            fw = open('gatk_snp_indel.sh','w').write(out2)
        elif combine == 'Y' and tool == 'samtools+gatk':
            fw = open('samtools_snp_indel.sh','w').write(out1)
            fw = open('gatk_snp_indel.sh','w').write(out2)
            fw = open('samtools_gatk_combine.sh','w').write(out3)

if __name__ == "main":
    examplelog = """EXAMPLES:
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g  
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -ct samtools -mode SNP+Indel
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -ct gatk4 -bqsr /root/my_data/known_sites/ -vqsr /root/my_data/resources/
    -pr 'Y' -s1 Tu4 -s2 Nm_35 -mode SNP+Indel
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -ct gatk4 -pr 'Y' -s1 Tu4 -s2 Nm_35 -mode SNP+Indel
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -ct control-freec -pr 'N' -mode CNV
    
    """
    parser = argparse.ArgumentParser(description='Variation pipline v1.0',
                                     epilog=examplelog,
                                     add_help=False)
    general = parser.add_argument_group(title='General options')
    general.add_argument('-h', '--help', action="help",
                         help="show the help and exit")
    general.add_argument('-i', '--input', type=str,
                         help="The input directory of clean reads")
    general.add_argument('-o', '--output', type=str,
                         help="The output directory")
    general.add_argument('-r', '--reference', type=str,
                         help="The fasta file of reference")
    general.add_argument('-g', '--gff3', type=str,
                         help="The gff3 file of reference for variation sites annotation")
    mapref = parser.add_argument_group(title='Mapping options')
    mapref.add_argument('-sp', '--seq_platform',type=str, default='ngs',choices=['ngs','tgs'],
                        help='Reads sequencing types, ngs for illumina short reads; tgs for Pacbio or ONT long reads')
    mapref.add_argument('-mt', '--maptool', type=str, default='BWA', choices=['BWA','Minimap2','NGML'],
                        help='Choose an alignment tool,illumina:BWA;Pacbio/ONT:Minimap2/NGML')
    mapref.add_argument('-mp', '--maptool_parameters',type=str,default='',
                        help="Set parameters for alignment tools")
    mutation = parser.add_argument_group(title='Detection variation (SNP/Indel/SV/CNV) options')

    mutation.add_argument('-ct','--calltool',type=str,default='samtools',choices=['samtools','gatk4','samtools+gakt4',
                                                                                  'breakdancer','crest',
                                                                                  'cnvnator','control-freec','sniffles'],
                          help='Choose a detection pipeline for SNP/Indel/SV/CNV')
    mutation.add_argument('-mode', '--mode', type=str, default='SNP+Indel',choices=['SNP+Indel','SV+CNV','SV','CNV'],
                          help='Mutation types')
    mutation.add_argument('-bqsr','--BQSR',type=str,default='',
                          help='If choose BQSR parameter(GATK4 only), you should provide a directory that only contain \
                          known sites')
    mutation.add_argument('-vqsr','--VQSR',type=str,default='',
                          help='Choose a method for mutations quality control(GATK4 only), if set this parameters, \
                               please provide a directory that only contain resource files, if not hard_filtering will \
                               be used')
    mutation.add_argument('-pr', '--pairs', type=str,default='N',choices=['N','Y'],
                          help='Detect germline and somatic mutation, default is N (Not to detect)')
    mutation.add_argument('-s1', '--sample',type=str,
                          help='The target(tumor/disease) sample name')
    mutation.add_argument('-s2', '--control', type=str,
                           help='The control(normal) sample name')




