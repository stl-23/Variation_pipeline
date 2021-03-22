import os
import re
from tools import parsering, getmyconfig


samtools = getmyconfig.getConfig('Variation', 'samtools')

class Mapping(object):
    def __init__(self, maptools, inputs, outputs, refs, parameters,platform,type):
        self.maptools = maptools
        self.inputs = inputs
        self.outputs = outputs
        self.refs = refs
        self.parameters = parameters
        self.platform = platform
        self.type = type
    def ngs(self):
        input_path = os.path.abspath(self.inputs) + '/'
        out_path = os.path.abspath(self.outputs) + '/'
        #lst = os.listdir(input_path)
        outfile = []
        samples, input1, input2, output1, output2 = parsering.parse_short_read_dir(input_path,out_path,self.type)
        if self.type == 'PE':
            for sample,index in enumerate(samples):
                if re.search('-R',self.parameters):
                    outfile.append('{bwa} {parameters} {ref} {input1} {input2} | {samtools} view -bhS - > {out_path}/{sample}.bam \
                                   && {samtools} sort -n -o {out_path}/{sample}.namesort.bam {out_path}/{sample}.bam \
                                   && {samtools} fixmate -m {out_path}/{sample}.namesort.bam {out_path}/{sample}.fixmate.bam \
                                   && {samtools} sort {out_path}/{sample}.fixmate.bam -o {out_path}/{sample}.pos.sort \
                                   && {samtools} markdup -r {out_path}/{sample}.pos.sort.bam {out_path}/{sample}.rmdup.bam \
                                   && {samtools} index {out_path}/{sample}.rmdup.bam'.format(
                        bwa=self.maptools,parameters=self.parameters,ref=self.refs,input1=input1,input2=input2,
                        samtools=samtools,out_path=out_path,sample=sample
                        ))
                else:
                    self.parameters = self.parameters+"-R '@RG\tID:{sample}\tPL:ILLUMINA\tSM:{sample}'".format(sample=sample)
                    outfile.append('{bwa} {parameters} {ref} {input1} {input2} | {samtools} view -bhS - > {out_path}/{sample}.bam \
                                   && {samtools} sort -n -o {out_path}/{sample}.namesort.bam {out_path}/{sample}.bam \
                                   && {samtools} fixmate -m {out_path}/{sample}.namesort.bam {out_path}/{sample}.fixmate.bam \
                                   && {samtools} sort {out_path}/{sample}.fixmate.bam -o {out_path}/{sample}.pos.sort \
                                   && {samtools} markdup -r {out_path}/{sample}.pos.sort.bam {out_path}/{sample}.rmdup.bam \
                                   && {samtools} index {out_path}/{sample}.rmdup.bam'.format(
                            bwa=self.maptools, parameters=self.parameters, ref=self.refs, input1=input1, input2=input2,
                            samtools=samtools, out_path=out_path, sample=sample
                        ))
        elif self.type == 'SE':
            for sample, index in enumerate(samples):
                if re.search('-R',self.parameters):
                    outfile.append('{bwa} {parameters} {ref} {input1} | {samtools} view -bhS - > {out_path}/{sample}.bam \
                    && {samtools} sort {out_path}/{sample}.bam -o {out_path}/{sample}.sorted \
                    && {samtools} index {out_path}/{sample}.sorted \
                    && {samtools} rmdup -s {out_path}/{sample}.sorted.bam {out_path}/{sample}.rmdup.bam'.format(
                        bwa=self.maptools, parameters=self.parameters, ref=self.refs, input1=input1,
                        samtools=samtools, out_path=out_path, sample=sample
                    ))
                else:
                    self.parameters = self.parameters + "-R '@RG\tID:{sample}\tPL:ILLUMINA\tSM:{sample}'".format(
                        sample=sample)
                    outfile.append('{bwa} {parameters} {ref} {input1} | {samtools} view -bhS - > {out_path}/{sample}.bam \
                    && {samtools} sort {out_path}/{sample}.bam -o {out_path}/{sample}.sorted \
                    && {samtools} index {out_path}/{sample}.sorted \
                    && {samtools} rmdup -s {out_path}/{sample}.sorted.bam {out_path}/{sample}.rmdup.bam'.format(
                        bwa=self.maptools, parameters=self.parameters, ref=self.refs, input1=input1,
                        samtools=samtools, out_path=out_path, sample=sample
                    ))
        return outfile

    def tgs_minimap2(self):
        input_path = os.path.abspath(self.inputs) + '/'
        out_path = os.path.abspath(self.outputs) + '/'
        #lst = os.listdir(input_path)
        outfile = []
        samples, input1 = parsering.parse_long_read_dir(input_path)
        for sample, index in enumerate(samples):
            outfile.append('{minimap2} {parameters} {ref} {input1} | {samtools} view -bhS - > {out_path}/{sample}.bam \
                           && {samtools} sort {out_path}/{sample}.bam -o {out_path}/{sample}.sorted \
                            && {samtools} index {out_path}/{sample}.sorted'.format(
                            minimap2=self.maptools, parameters=self.parameters, ref=self.refs, input1=input1,
                        samtools=samtools, out_path=out_path, sample=sample
            ))
        return outfile

    def tgs_ngml(self):
        input_path = os.path.abspath(self.inputs) + '/'
        out_path = os.path.abspath(self.outputs) + '/'
        #lst = os.listdir(input_path)
        outfile = []
        samples, input1 = parsering.parse_long_read_dir(input_path)
        for sample, index in enumerate(samples):
            outfile.append('{ngml} {parameters} {ref} {input1} | {samtools} view -bhS - > {out_path}/{sample}.bam \
                                       && {samtools} sort {out_path}/{sample}.bam -o {out_path}/{sample}.sorted \
                                        && {samtools} index {out_path}/{sample}.sorted'.format(
                ngml=self.maptools, parameters=self.parameters, ref=self.refs, input1=input1,
                samtools=samtools, out_path=out_path, sample=sample
            ))
        return outfile
