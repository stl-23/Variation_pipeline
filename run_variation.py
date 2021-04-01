#!/usr/bin/env python
import os
import sys
import argparse
import mapping
import annotation
import ngs_vars
import tgs_vars
from vartools import parsering

#script_path = os.path.abspath("./")
#make_index = os.path.join(script_path, '/tools/index.sh')

def rw_mapping():
    try:
        input_path = os.path.abspath(inputs_dir) + '/'
        out_path = os.path.abspath(outputs_dir) + '/mapping_results'
        os.system('mkdir -p {0}'.format(out_path))
        # lst = os.listdir(input_path)
        outfile = []
        global samples
        samples = set()
        if platform == "ngs":
            samples = parsering.parse_short_read_dir(input_path, out_path, seqtype)[0]

            outfile = mapping.Ngs(maptool, input_path, out_path,ref,maptool_parameters,seqtype).ngs()
            for i,s in enumerate(outfile):
                fw = open('ngs'+i+'.bwa.sh','w').write(s)

        elif platform == 'tgs':
            samples = parsering.parse_long_read_dir(input_path)[0]
            if maptool == 'minimap2':
                outfile = mapping.Tgs(maptool,input_path,out_path,ref,maptool_parameters).tgs_minimap2()
                for i,s in enumerate(outfile):
                    fw = open('tgs'+i+'.minimap2.sh','w').write(s)
            elif maptool == 'ngml':
                outfile = mapping.Tgs(maptool,input_path,out_path,ref,maptool_parameters).tgs_ngml()
                for i,s in enumerate(outfile):
                    fw = open('tgs'+i+'.ngml.sh','w').write(s)

    except Exception as e:
        print('Error files or mapping tools')

def rw_call_var():
    try:
        input_path = os.path.abspath(outputs_dir) + '/mapping_results'
        out_path = os.path.abspath(outputs_dir) + '/var_results'
        os.system('mkdir -p {0}'.format(out_path))
        # lst = os.listdir(input_path)
        input = [input_path + x for x in list(samples)]
        samples = [out_path + x for x in list(samples)]
        outfile = []
        if platform == 'ngs':
            if mode == 'SNP_Indel':
                if calltool == 'samtools':
                    if v_calling == 'single':
                        for sample in samples:
                            outfile.append(ngs_vars.snp_indel_samtools(ref, input, sample))
                elif calltool == 'gatk4':
                    if bqsr == '' and vqsr == '':
                        for sample in samples:
                            outfile.append(ngs_vars.snp_indel_gatk(ref, input, sample, bqsr, ))
                elif calltool == 'samtools+gatk4':
            elif mode == 'SV':
            elif mode == 'CNV':
            elif mode == 'SNP_Indel_Somatic':
            elif mode == 'SV_Somatic':
            elif mode == 'CNV_Somatic':

        elif platform == 'tgs':
            if mode == 'SNP_Indel':
            elif mode == 'SV':

def rw_annotation():
    try:
        input_path = os.path.abspath(outputs_dir) + '/var_results'
        out_path = os.path.abspath(outputs_dir) + '/annotation_results'
        os.system('mkdir -p {0}'.format(out_path))


if __name__ == '__main__':
    examplelog = """EXAMPLES:
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g  
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -ct samtools -mode SNP+Indel
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -ct gatk4 -bqsr /root/my_data/known_sites/ -vqsr /root/my_data/resources/
    -pr 'Y' -s1 Tu4 -s2 Nm_35 -mode SNP+Indel
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -ct gatk4 -pr 'Y' -s1 Tu4 -s2 Nm_35 -mode SNP+Indel
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -ct control-freec -pr 'N' -mode CNV -ot WGS 
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp tgs -mt minimap2 -ct sniffles -mode SV 
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
                        help='Reads sequencing platform, ngs for illumina short reads; tgs for Pacbio or ONT long reads')
    mapref.add_argument('-st', '--seq_type', type=str, default='PE', choices=['PE', 'SE'],
                        help='Sequencing type, paired-end/mate-pair, or single end')
    mapref.add_argument('-mt', '--maptool', type=str, default='BWA', choices=['BWA','Minimap2','NGML'],
                        help='Choose an alignment tool,illumina:BWA;Pacbio/ONT:Minimap2/NGML')
    mapref.add_argument('-mp', '--maptool_parameters',type=str,default='',
                        help="Set parameters for alignment tools")
    mutation = parser.add_argument_group(title='Detection variation (SNP/Indel/SV/CNV) options')

    mutation.add_argument('-ct','--calltool',type=str,default='samtools',choices=['samtools','gatk4','samtools+gakt4',
                                                                                  'breakdancer','crest',
                                                                                  'cnvnator','control-freec','sniffles'],
                          help='Choose a detection pipeline for SNP/Indel/SV/CNV')
    mutation.add_argument('-mode', '--mode', type=str, default='SNP_Indel',choices=['SNP_Indel','SV','CNV',
                                                                              'Germline','Somatic'],
                          help='Mutation types')
    mutation.add_argument('-bqsr','--BQSR',type=str,default='',
                          help='If choose BQSR parameter(GATK4 only), you should provide a directory that only contain \
                          known sites')
    mutation.add_argument('-vqsr','--VQSR',type=str,default='',
                          help='Choose a method for mutations quality control(GATK4 only), if set this parameters, \
                               please provide a directory that only contain resource files, if not hard_filtering will \
                               be used')
    mutation.add_argument('-vc','--variation_calling', type=str,default='single',choices=['single','join'],
                           help='Calling a group of samples together(join calling) or Variant calling with a single \
                           sample only(single sample calling), default is single')
    mutation.add_argument('-pr', '--pairs', type=str,default='N',choices=['N','Y'],
                          help='Detect germline and somatic mutation, default is N (Not to detect)')
    mutation.add_argument('-s1', '--sample',type=str,
                          help='The target(tumor/disease) sample name')
    mutation.add_argument('-s2', '--control', type=str,
                           help='The control(normal) sample name')
    mutation.add_argument('-sg', '--strategy', type=str, default='WGS',choices=['WGS','WES'],
                          help='WGS or WES')
    args = parser.parse_args()
    inputs_dir = args.input
    outputs_dir = args.outputs
    ref = args.reference
    gff3 = args.gff3
    platform = args.seq_platform
    seqtype = args.seq_type
    maptool = args.maptool
    maptool_parameters = args.maptool_parameters
    calltool = args.calltool
    mode = args.mode
    bqsr = args.BQSR
    vqsr = args.VQSR
    v_calling = args.variation_calling
    pairs = args.pairs
    sample = args.sample
    control = args.control
    strategy = args.strategy

    rw_mapping()
    rw_call_var()
    rw_annotation()



