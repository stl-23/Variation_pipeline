import os
import re
import argparse
from varition import getmyconfig
from parsedic import parsering


class Mapping(object):
    def __init__(self, maptools, inputs, outputs, refs, parameters,platform,type):
        self.maptools = maptools
        self.inputs = inputs
        self.outputs = outputs
        self.refs = refs
        self.parameters = parameters
        self.platform = platform
        self.type = type
    def parse(self):
        input_path = os.path.abspath(self.inputs) + '/'
        out_path = os.path.abspath(self.outputs) + '/'
        lst = os.listdir(input_path)
        outfile = []
        if self.platform == 'ngs':
            samples, input1, input2, output1, output2 = parsering.parse_short_read_dir(input_path,out_path,self.type)
            if self.type == 'PE':
                for sample,index in enumerate(samples):
                    if re.search('-R',self.parameters):
                        outfile.append('{bwa} {parameters} {ref} {input1} {input2} | {samtools}')
                    else:
                        self.parameters = self.parameters+"-R '@RG\tID:{sample}\tPL:Illumina\tSM:{sample}'".format(sample=sample)


        elif self.platform == 'tgs':
            tgs_seq_suffix = ['.fa', '.fasta', '.fastq', '.fq', '.fq.gz', '.fastq.gz']
            tgs_lst = [file for file in lst for suffix in tgs_seq_suffix if file.endswith(suffix)]
            samples = [i.replace('.fa','').replace('.fasta','').replace('.fastq','').replace('.fq','').replace('.fq.gz','').replace('.fastq.gz','') for i in tgs_lst]
            samples = sorted(samples)

    def ngs(self):
        parse()
    def tgs(self):
