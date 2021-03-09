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
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq.gz')
                input2.append(input_path + sample + '_2.fq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')

        elif lst[0].endswith('.fq'):
            samples = [i.replace('_1.fq', '').replace('_2.fq', '') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq')
                input2.append(input_path + sample + '_2.fq')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')

        elif lst[0].endwith('.fastq.gz'):
            samples = [i.replace('_1.fastq.gz', '').replace('_2.fastq.gz', '') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq.gz')
                input2.append(input_path + sample + '_2.fastq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')

        elif lst[0].endwith('.fastq'):
            samples = [i.replace('_1.fastq', '').replace('_2.fastq', '') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq')
                input2.append(input_path + sample + '_2.fastq')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')

    elif seq_type == 'SE':
        if lst[0].endswith('.fq.gz'):
            samples = [i.replace('.fq.gz') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')

        elif lst[0].endswith('.fq'):
            samples = [i.replace('.fq') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq')
                output1.append(out_path + sample + '_1.clean.fq.gz')

        elif lst[0].endwith('.fastq.gz'):
            samples = [i.replace('.fastq.gz') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')

        elif lst[0].endwith('.fastq'):
            samples = [i.replace('.fastq') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq')
                output1.append(out_path + sample + '_1.clean.fq.gz')


    return samples, input1, input2, output1, output2

def parse_pacbio_read_dir(inputs, outs, mt_type='Sequel'):
    input_path = os.path.abspath(inputs) + '/'
    out_path = os.path.abspath(outs) + '/'
    lst = os.listdir(input_path)
    input1=[]
    input2=[]
    input3=[]
    output1=[]
    samples=[]
    try:
        seq_suffix = ['.subreads.bam','1.bax.h5','.2.bax.h5','.3.bax.h5']
        lst = [file for file in lst for suffix in seq_suffix if file.endswith(suffix)]
    except Exception:
        print("No such file or directory or wrong file format,must be .subreads.bam/.1/2/3.bax.h5")
        exit()
    if mt_type == 'Sequel':
        samples = [i.replace('.subreads.bam') for i in lst if i.endswith('.subreads.bam')]  # or .subreadset.xml
        samples = sorted(samples)
        for sample in samples:
            input1.append(input_path + sample)
            output1.append(out_path + sample)
    elif mt_type == 'RS':
        samples = [i.replace('.1.bax.h5','').replace('.2.bax.h5','').replace('.3.bax.h5','') for i in lst if i.endswith('.bax.h5')]  ## one bas.h5 file and three bax.h5 files in each sample run
        samples = sorted(samples)
        for sample in samples:
            input1.append(input_path + sample + '.1.bax.h5')
            input2.append(input_path + sample + '.2.bax.h5')
            input3.append(input_path + sample + '.3.bax.h5')
            output1.append(out_path + sample)

    return samples,inuput1,input2,input3,output1

def parse_nanopore_read_dir(inputs, outs):
    input_path = os.path.abspath(inputs) + '/'
    out_path = os.path.abspath(outs) + '/'
    lst = os.listdir(input_path)
    samples = [i.replace('.fastq.gz','') for i in lst if i.endswith('.fastq.gz')]
    samples = sorted(samples)
    input_samples = []
    for sample in samples:
        input_sample.append(input_path + sample + '.fastq.gz')
    return samples, input_samples
