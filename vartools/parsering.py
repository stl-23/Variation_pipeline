
import os

def parse_short_read_dir(inputs, outs, seq_type='PE'):
    input_path = os.path.abspath(inputs)
    out_path = os.path.abspath(outs)
    lst = os.listdir(input_path)
    input1 = []
    input2 = []
    output1 = []
    output2 = []
    seq_suffix = [seq + zz for seq in ['.fq', '.fastq'] for zz in ['', '.gz']]
    lst = [file for file in lst for s in seq_suffix if file.endswith(s)]
    if lst:
        samples = [ i.split('_1')[0] for i in lst if '_1' in i]
        suffix = [i.split('_1')[1] for i in lst if '_1' in i]
    else:
        raise IOError('No such file or directory or wrong file suffix (e.q. sample_1.fq.gz/_1.fastq.gz/_1.fq/_1.fastq)')
    if samples:
        if seq_type == 'PE':
            for index, sample in enumerate(samples):
                input1.append(os.path.join(input_path,sample+'_1'+suffix[index]))
                input2.append(os.path.join(input_path,sample+'_2'+suffix[index]))
                output1.append(os.path.join(out_path,sample+'_1'+suffix[index]))
                output2.append(os.path.join(out_path,sample+'_2'+suffix[index]))

        elif seq_type == 'SE':
            for index, sample in enumerate(samples):
                input1.append(os.path.join(input_path,sample+'_1'+suffix[index]))
                output1.append(os.path.join(out_path, sample + '_1' + suffix[index]))

    else:
        raise IOError('Wrong file name (e.q. sample_1.fq.gz/_1.fastq.gz/_1.fq/_1.fastq)')


    return samples, input1, input2, output1, output2

def parse_long_read_dir(inputs):  ## clean data
    input_path = os.path.abspath(inputs) + '/'
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

