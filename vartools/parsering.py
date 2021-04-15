
import os

def parse_short_read_dir(inputs, outs, seq_type='PE'):
    input_path = os.path.abspath(inputs) + '/'
    out_path = os.path.abspath(outs) + '/'
    lst = os.listdir(input_path)
    samples = []
    input1 = []
    input2 = []
    output1 = []
    output2 = []

    try:
        seq_suffix = ['.fq','.fq.gz','.fastq','.fastq.gz']
        lst = [file for file in lst for suffix in seq_suffix if file.endswith(suffix)]
    except Exception:
        print("No such file or directory or wrong file format,must be .fq/.fq.gz/.fastq/.fastq.gz")
        exit()
    if seq_type == 'PE':
        if lst[0].endswith('.fq.gz'):
            samples = [i.replace('_1.fq.gz', '').replace('_2.fq.gz', '') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq.gz')
                input2.append(input_path + sample + '_2.fq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')

        elif lst[0].endswith('.fq'):
            samples = [i.replace('_1.fq', '').replace('_2.fq', '') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq')
                input2.append(input_path + sample + '_2.fq')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')

        elif lst[0].endswith('.fastq.gz'):
            samples = [i.replace('_1.fastq.gz', '').replace('_2.fastq.gz', '') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq.gz')
                input2.append(input_path + sample + '_2.fastq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')

        elif lst[0].endswith('.fastq'):
            samples = [i.replace('_1.fastq', '').replace('_2.fastq', '') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq')
                input2.append(input_path + sample + '_2.fastq')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')

    elif seq_type == 'SE':
        if lst[0].endswith('.fq.gz'):
            samples = [i.replace('.fq.gz') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')

        elif lst[0].endswith('.fq'):
            samples = [i.replace('.fq') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq')
                output1.append(out_path + sample + '_1.clean.fq.gz')

        elif lst[0].endswith('.fastq.gz'):
            samples = [i.replace('.fastq.gz') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')

        elif lst[0].endswith('.fastq'):
            samples = [i.replace('.fastq') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq')
                output1.append(out_path + sample + '_1.clean.fq.gz')


    return samples, input1, input2, output1, output2

def parse_long_read_dir(inputs):  ## clean data
    input_path = os.path.abspath(inputs) + '/'
    #out_path = os.path.abspath(outs) + '/'
    lst = os.listdir(input_path)
    input1 = []
    tgs_seq_suffix = ['.fa', '.fasta', '.fastq', '.fq', '.fa.gz', 'fasta.gz', '.fq.gz', '.fastq.gz']
    tgs_lst = [file for file in lst for suffix in tgs_seq_suffix if file.endswith(suffix)]
    samples = [i.replace('.fastq.gz', '').replace('.fasta.gz','').replace('.fastq', '').replace('.fasta', '').replace('fa.gz','').replace('.fq.gz','').replace('.fa', '').replace('.fq', '') for i in tgs_lst]
    samples = set(samples)
    if tgs_lst[0].endswith('.fa'):
        for sample in samples:
            input1.append(input_path + sample + '.fa')
    elif tgs_lst[0].endswith('.fasta'):
        for sample in samples:
            input1.append(input_path + sample + '.fasta')
    elif tgs_lst[0].endswith('.fastq'):
        for sample in samples:
            input1.append(input_path + sample + '.fastq')
    elif tgs_lst[0].endswith('.fq'):
        for sample in samples:
            input1.append(input_path + sample + '.fq')
    elif tgs_lst[0].endswith('.fa.gz'):
        for sample in samples:
            input1.append(input_path + sample + '.fa.gz')
    elif tgs_lst[0].endswith('.fasta.gz'):
        for sample in samples:
            input1.append(input_path + sample + '.fasta.gz')
    elif tgs_lst[0].endswith('.fq.gz'):
        for sample in samples:
            input1.append(input_path + sample + '.fq.gz')
    elif tgs_lst[0].endswith('.fastq.gz'):
        for sample in samples:
            input1.append(input_path + sample + '.fastq.gz')

    return samples, input1
'''
def common_parse(inputs,outs,*suffix):
    input_path = os.path.abspath(inputs) + '/'
    out_path = os.path.abspath(outs) + '/'
    lst = os.listdir(input_path)
    input1 = []
    files = [file for file in lst for s in suffix if file.endswith(s)]
    samples = set([f.replace(x) for f in files for s in suffix])
    
'''
