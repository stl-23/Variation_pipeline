import os
import re
from vartools import parsering, getmyconfig


samtools = getmyconfig.getConfig('Variation', 'samtools')
#picard = getmyconfig.getConfig(('Variation','picard'))
gatk4 = getmyconfig.getConfig('Variation','gatk4')

class Mapping(object):
    def __init__(self, maptools, inputs, outputs, refs, parameters):
        self.maptools = maptools
        self.inputs = inputs
        self.outputs = outputs
        self.refs = refs
        self.parameters = parameters
    def parse(self):
        input_path = os.path.abspath(self.inputs) + '/'
        out_path = os.path.abspath(self.outputs) + '/'
        #lst = os.listdir(input_path)
        outfile = []
        return self.maptools,self.parameters,self.refs,input_path,out_path,outfile

class Ngs(Mapping):
    def __init__(self,maptools, inputs, outputs, refs, parameters,type):
        super(Ngs,self).__init__(maptools, inputs, outputs, refs, parameters)
        self.type = type
        self.maptools,self.parameters,self.refs,self.input_path,self.out_path,self.outfile = super().parse()
    def ngs(self):
        samples, input1, input2, output1, output2 = parsering.parse_short_read_dir(self.input_path,self.out_path,self.type)
        if self.type == 'PE':
            for index,sample in enumerate(samples):
                if re.search('-R',self.parameters):  ## bwa mem -R '@RG\tID:foo\tPL:ILLUMINA\tSM:foo' ref.fa read1.fq read2.fq
                    #self.outfile.append('{bwa} {parameters} {ref} {input1} {input2} | {samtools} view -bhS - > {out_path}/{sample}.bam \
                    #    && {samtools} sort -n -o {out_path}/{sample}.namesort.bam {out_path}/{sample}.bam \
                    #    && {samtools} fixmate -m {out_path}/{sample}.namesort.bam {out_path}/{sample}.fixmate.bam \
                    #    && {samtools} sort {out_path}/{sample}.fixmate.bam -o {out_path}/{sample}.pos.sort \
                    #    && {samtools} markdup -r {out_path}/{sample}.pos.sort.bam {out_path}/{sample}.rmdup.bam \
                    #    && {samtools} index {out_path}/{sample}.rmdup.bam'.format(
                    #    bwa=self.maptools,parameters=self.parameters,ref=self.refs,input1=input1,input2=input2,
                    #    samtools=samtools,out_path=self.out_path,sample=sample
                    #    ))
                    self.outfile.append("""{bwa} {parameters} {ref} {input1} {input2} | {samtools} view -bhS - > {out_path}/{sample}.bam 
&& {samtools} sort {out_path}/{sample}.bam -o {out_path}/{sample}.pos.sort 
&& {gatk4} --java-options "-Xmx8G" MarkDuplicates -I {sample}.pos.sort.bam -O {sample}.rmdup.bam -M {sample}.rmdup.txt --REMOVE_DUPLICATES true 
&& {gatk4} --java-options "-Xmx8G" BuildBamIndex -I {sample}.rmdup.bam
&& {samtools} flagstat {sample}.rmdup.bam > {sample}.rmdup.stat""".format(
                        bwa=self.maptools,parameters=self.parameters,ref=self.refs,input1=input1[index],input2=input2[index],
                        samtools=samtools,out_path=self.out_path,sample=sample,gatk4=gatk4,
                        ))
                else:
                    self.parameters = self.parameters+" -R '@RG\tID:{sample}\tLB:{sample}\tPL:ILLUMINA\tSM:{sample}'".format(sample=sample)
                    #self.outfile.append('{bwa} {parameters} {ref} {input1} {input2} | {samtools} view -bhS - > {out_path}/{sample}.bam \
                    #            && {samtools} sort -n -o {out_path}/{sample}.namesort.bam {out_path}/{sample}.bam \
                    #           && {samtools} fixmate -m {out_path}/{sample}.namesort.bam {out_path}/{sample}.fixmate.bam \
                    #            && {samtools} sort {out_path}/{sample}.fixmate.bam -o {out_path}/{sample}.pos.sort \
                    #            && {samtools} markdup -r {out_path}/{sample}.pos.sort.bam {out_path}/{sample}.rmdup.bam \
                    #            && {samtools} index {out_path}/{sample}.rmdup.bam'.format(
                    #    bwa=self.maptools, parameters=self.parameters, ref=self.refs, input1=input1, input2=input2,
                    #    samtools=samtools, out_path=self.out_path, sample=sample
                    #    ))
                    self.outfile.append("""{bwa} {parameters} {ref} {input1} {input2} | {samtools} view -bhS - > {out_path}/{sample}.bam 
&& {samtools} sort {out_path}/{sample}.bam -o {out_path}/{sample}.pos.sort 
&& {gatk4} --java-options "-Xmx8G" MarkDuplicates -I {sample}.pos.sort.bam -O {sample}.rmdup.bam -M {sample}.rmdup.txt --REMOVE_DUPLICATES true 
&& {gatk4} --java-options "-Xmx8G" BuildBamIndex -I {sample}.rmdup.bam 
&& {samtools} flagstat {sample}.rmdup.bam > {sample}.rmdup.stat""".format(
                        bwa=self.maptools,parameters=self.parameters,ref=self.refs,input1=input1[index],input2=input2[index],
                        samtools=samtools,out_path=self.out_path,sample=sample,gatk4=gatk4
                        ))
        elif self.type == 'SE':
            for index, sample in enumerate(samples):
                if re.search('-R',self.parameters):
                    self.outfile.append("""{bwa} {parameters} {ref} {input1} | {samtools} view -bhS - > {out_path}/{sample}.bam 
&& {samtools} sort {out_path}/{sample}.bam -o {out_path}/{sample}.sorted 
&& {samtools} index {out_path}/{sample}.sorted
&& {samtools} rmdup -s {out_path}/{sample}.sorted.bam {out_path}/{sample}.rmdup.bam
&& {samtools} flagstat {sample}.rmdup.bam > {sample}.rmdup.stat""".format(
                        bwa=self.maptools, parameters=self.parameters, ref=self.refs, input1=input1[index],
                        samtools=samtools, out_path=self.out_path, sample=sample
                        ))
                else:
                    self.parameters = self.parameters + " -R '@RG\tID:{sample}\tLB:{sample}\tPL:ILLUMINA\tSM:{sample}'".format(
                            sample=sample)
                    self.outfile.append("""{bwa} {parameters} {ref} {input1} | {samtools} view -bhS - > {out_path}/{sample}.bam
&& {samtools} sort {out_path}/{sample}.bam -o {out_path}/{sample}.sorted
&& {samtools} index {out_path}/{sample}.sorted
&& {samtools} rmdup -s {out_path}/{sample}.sorted.bam {out_path}/{sample}.rmdup.bam
&& {samtools} flagstat {sample}.rmdup.bam > {sample}.rmdup.stat""".format(
                        bwa=self.maptools, parameters=self.parameters, ref=self.refs, input1=input1[index],
                        samtools=samtools, out_path=self.out_path, sample=sample
                        ))
        return self.outfile
class Tgs(Mapping):
    def __init__(self,maptools, inputs, outputs, refs, parameters):
        super(Tgs,self).__init__( maptools, inputs, outputs, refs, parameters)
        self.maptools,self.parameters,self.refs,self.input_path,self.out_path,self.outfile = super().parse()
    def tgs_minimap2(self):
        samples, input1 = parsering.parse_long_read_dir(self.input_path)
        for index, sample in enumerate(samples):
            self.outfile.append("""{minimap2} {parameters} {ref} {input1} | {samtools} view -bhS - > {out_path}/{sample}.bam 
&& {samtools} sort {out_path}/{sample}.bam -o {out_path}/{sample}.sorted
&& {samtools} index {out_path}/{sample}.sorted
&& {samtools} flagstat {sample}.sorted.bam > {sample}.sorted.stat""".format(
                            minimap2=self.maptools, parameters=self.parameters, ref=self.refs, input1=input1[index],
                        samtools=samtools, out_path=self.out_path, sample=sample
            ))
        return self.outfile

    def tgs_ngml(self):
        samples, input1 = parsering.parse_long_read_dir(self.input_path)
        for index, sample in enumerate(samples):
            self.outfile.append("""{ngml} {parameters} {ref} {input1} | {samtools} view -bhS - > {out_path}/{sample}.bam
&& {samtools} sort {out_path}/{sample}.bam -o {out_path}/{sample}.sorted
&& {samtools} index {out_path}/{sample}.sorted
&& {samtools} flagstat {sample}.sorted.bam > {sample}.sorted.stat""".format(
                ngml=self.maptools, parameters=self.parameters, ref=self.refs, input1=input1[index],
                samtools=samtools, out_path=self.out_path, sample=sample
            ))
        return self.outfile
