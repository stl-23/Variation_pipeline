#!/usr/bin/env python
import os
import sys
import argparse
import mapping
import annotation
import ngs_vars
import tgs_vars
from vartools import parsering, merge_vcf_gcvf,variation_qc,somatic_detection

#script_path = os.path.abspath("./")
#make_index = os.path.join(script_path, '/tools/index.sh')

def write_mapping():
    try:
        input_path = os.path.abspath(inputs_dir) + '/'
        out_path = os.path.join(os.path.abspath(outputs_dir), '/mapping_results/')
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


def write_call_var():
    input_path = os.path.join(os.path.abspath(outputs_dir), 'mapping_results/')
    out_path = os.path.join(os.path.abspath(outputs_dir), 'var_results/')
    if buildver == 'hg19':
        genomicsdb = os.path.join(os.path.abspath(inputs_dir), 'genomicsdb/hg19/')
        chr_list = os.path.join(os.path.abspath(inputs_dir), '/hg19_chr.list')
    elif buildver == 'hg38':
        genomicsdb = os.path.join(os.path.abspath(inputs_dir), 'genomicsdb/hg38/')
        chr_list = os.path.join(os.path.abspath(inputs_dir), 'hg38_chr.list')
    else:
        genomicsdb = os.path.join(os.path.abspath(inputs_dir), 'genomicsdb/{0}/'.format(buildver))
        chr_list = os.path.join(os.path.abspath(inputs_dir), '{0}_chr.list'.format(buildver))

    tmp_dir = os.path.join(os.path.abspath(outputs_dir), 'tmp_dir/')
    os.system('mkdir -p {0} {1}'.format(out_path,tmp))
    # lst = os.listdir(input_path)
    path_input = [input_path + x for x in list(samples)]
    path_output = [out_path + x for x in list(samples)]
    for s in samples:
        fw = open('map_file.list','w').write(s+'\t'+s+'.g.vcf')
    if platform == 'ngs':
        if mode == 'SNP_Indel':
            if callpipe == 'samtools':
                cmd_call = []
                cmd_merge = ''
                if v_calling == 'single': ##  single sample calling
                    for index,sample in enumerate(path_output):
                        cmd_call.append(ngs_vars.snp_indel_samtools(ref, path_input[index], sample))
                    for index, per_cmd in enumerate(cmd_call):
                        fw = open(samples[index]+'ngs_samtools_call.sh','w').write(per_cmd)
                elif v_calling == 'join':
                    for index,sample in enumerate(path_output): ## join calling
                        cmd_call.append(ngs_vars.snp_indel_samtools(ref, path_input[index], sample))
                    for index, per_cmd in enumerate(cmd_call):
                        fw = open(samples[index] + 'ngs_samtools_call.sh', 'w').write(per_cmd)
                    files_snp = [i+'.samtools.raw.SNP.vcf' for i in path_output]
                    files_indel = [i+'.samtools.raw.INDEL.vcf' for i in path_output]
                    cmd_merge = merge_vcf_gcvf.merge(files_snp,'vcf',out_path+'all_samples_snp',[])
                    cmd_merge += merge_vcf_gcvf.merge(files_indel,'vcf',out_path+'all_samples_indel',[])
                    fw = open('ngs_samtools_merge.sh','w').write(cmd_merge)
            elif callpipe == 'gatk4':
                if bqsr_dir and vqsr_dir: ## BQSR and VQSR
                    cmd_call = []
                    cmd_merge = ''
                    if v_calling == 'single': ##  single sample calling
                        cmd_vqsr = []
                        for index,sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'T','F',bqsr_dir))
                            cmd_vqsr.append(variation_qc.vqsr(ref,path_output[index]+'.vcf',vqsr_dir,path_output[index]))     ### VQSR ###
                        for i, s in enumerate(samples):
                            fw = open(s+'ngs_gatk_bqsr_call_vqsr.sh','w').write(cmd_call[i]+'\n'+cmd_vqsr[i]+'\n')
                    elif v_calling == 'join':  ## join calling
                        cmd_vqsr = ''
                        for sample in path_output:
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample, 'T','T',bqsr_dir))
                        for index, per_cmd in enumerate(cmd_call):
                            fw = open(samples[index] + 'ngs_gatk_bqsr_call.sh', 'w').write(per_cmd)
                        files_g_vcf = [i+'.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gcvf.merge(files_g_vcf,'gvcf',out_path+'all_samples',
                                                            [4,ref,genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir])
                        cmd_vqsr = variation_qc.vqsr(ref,out_path+'all_samples.combined.vcf',vqsr_dir,out_path+"all_samples")
                        fw = open('ngs_gatk_bqsr_merge_vqsr.sh','w').write(cmd_merge+'\n'+cmd_vqsr)

                elif bqsr_dir and not vqsr_dir: ## BQSR and Hard filtering
                    cmd_call = []
                    cmd_merge = ''
                    if v_calling == 'single':  ##  single sample calling
                        cmd_hard = []
                        for index, sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'T', 'F', bqsr_dir))
                            cmd_hard.append(variation_qc.hard_filter(ref, path_output[index]+'.vcf', path_output[index]))
                        for i, s in enumerate(samples):
                            fw = open(s + 'ngs_gatk_bqsr_call_hard_filter.sh', 'w').write(cmd_call[i] + '\n' + cmd_hard[i] + '\n')
                    elif v_calling == 'join':  ## join calling
                        cmd_hard = ''
                        for sample in path_output:
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample, 'T', 'T', bqsr_dir))
                        for index, per_cmd in enumerate(cmd_call):
                            fw = open(samples[index] + 'ngs_gatk_bqsr_call.sh', 'w').write(per_cmd)
                        files_g_vcf = [i + '.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gcvf.merge(files_g_vcf, 'gvcf', out_path + 'all_samples',
                                                    [4, ref, genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir])
                        cmd_hard = variation_qc.hard_filter(ref, out_path + 'all_samples.combined.vcf', out_path + "all_samples")
                        fw = open('ngs_gatk_bqsr_merge_hard_filter.sh', 'w').write(cmd_merge + '\n' + cmd_hard)

                elif not bqsr_dir and not vqsr_dir: ## Hard filtering only
                    cmd_call = []
                    cmd_merge = ''
                    if v_calling == 'single':  ##  single sample calling
                        cmd_hard = []
                        for index, sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'F', 'F', bqsr_dir))
                            cmd_hard.append(variation_qc.hard_filter(ref, path_output[index] + '.vcf', path_output[index]))
                        for i, s in enumerate(samples):
                            fw = open(s + 'ngs_gatk_call_hard_filter.sh', 'w').write(cmd_call[i] + '\n' + cmd_hard[i] + '\n')
                    elif v_calling == 'join':  ## join calling
                        cmd_hard = ''
                        for sample in path_output:
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample, 'F', 'T', bqsr_dir))
                        for index, per_cmd in enumerate(cmd_call):
                            fw = open(samples[index] + 'ngs_gatk_call.sh', 'w').write(per_cmd)
                        files_g_vcf = [i + '.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gcvf.merge(files_g_vcf, 'gvcf', out_path + 'all_samples',
                                                         [4, ref, genomicsdb, chr_list, 1, 'map_file.list', 1, 26,
                                                          tmp_dir])
                        cmd_hard = variation_qc.hard_filter(ref, out_path + 'all_samples.combined.vcf',
                                                            out_path + "all_samples")
                        fw = open('ngs_gatk_merge_hard_filter.sh', 'w').write(cmd_merge + '\n' + cmd_hard)

            elif callpipe == 'samtools+gatk4':
                if v_calling == 'single':  ##  single sample calling
                    cmd_samtools_call = []
                    cmd_gatk4_call = []
                    cmd_combine = []
                    cmd_vqsr = []
                    cmd_hard = []
                    ###samtools pipeline###
                    for index,sample in enumerate(path_output):
                        cmd_samtools_call.append(ngs_vars.snp_indel_samtools(ref, path_input[index], sample))
                    for index, per_cmd in enumerate(cmd_samtools_call):
                        fw = open(samples[index]+'ngs_samtools_call.sh','w').write(per_cmd)
                    ###gatk4 pipeline###
                    if bqsr_dir and vqsr_dir:  ## BQSR and VQSR
                        for index, sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'T','F',bqsr_dir))
                            cmd_vqsr.append(variation_qc.vqsr(ref,path_output[index]+'.vcf',vqsr_dir,path_output[index]))
                        for i, s in enumerate(samples):
                            fw = open(s + 'ngs_gatk_bqsr_call_vqsr.sh', 'w').write(cmd_call[i] + '\n' + cmd_vqsr[i] + '\n')
                    elif bqsr_dir and not vqsr_dir: ## BQSR and Hard filtering
                        for index, sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'T', 'F', bqsr_dir))
                            cmd_hard.append(variation_qc.hard_filter(ref, path_output[index] + '.vcf', path_output[index]))
                        for i, s in enumerate(samples):
                            fw = open(s + 'ngs_gatk_call_hard_filter.sh', 'w').write(cmd_call[i] + '\n' + cmd_hard[i] + '\n')
                    elif not bqsr_dir and not vqsr_dir: ## Hard filtering only
                        for index, sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'F', 'F', bqsr_dir))
                            cmd_hard.append(variation_qc.hard_filter(ref, path_output[index] + '.vcf', path_output[index]))
                        for i, s in enumerate(samples):
                            fw = open(s + 'ngs_gatk_call_hard_filter.sh', 'w').write(cmd_call[i] + '\n' + cmd_hard[i] + '\n')
                    ###combine pipeline###
                    for index, sample in enumerate(path_output):
                        cmd_combine.append(ngs_vars.samtool_gatk_combine(path_output[index]))
                        for i, s in enumerate(samples):
                            fw = open(s + 'ngs_samtools_gatk_bqsr_vqsr_combine.sh', 'w').write(cmd_combine[i] + '\n')

                elif v_calling == 'join':  ## join calling
                    cmd_samtools_call = []
                    cmd_gatk4_call = []
                    cmd_merge = ''
                    cmd_combine = ''
                    cmd_vqsr = ''
                    cmd_hard = ''
                    ###samtools pipeline###
                    for index,sample in enumerate(path_output): ## join calling
                        cmd_samtools_call.append(ngs_vars.snp_indel_samtools(ref, path_inpu[index], sample))
                    for index, per_cmd in enumerate(cmd_samtools_call):
                        fw = open(samples[index] + 'ngs_samtools_call.sh', 'w').write(per_cmd)
                    files_snp = [i+'.samtools.raw.SNP.vcf' for i in path_output]
                    files_indel = [i+'.samtools.raw.INDEL.vcf' for i in path_output]
                    cmd_merge = merge_vcf_gcvf.merge(files_snp,'vcf',out_path+'all_samples_snp',[])
                    cmd_merge += merge_vcf_gcvf.merge(files_indel,'vcf',out_path+'all_samples_indel',[])
                    fw = open('ngs_samtools_merge.sh','w').write(cmd_merge)
                    ###gatk4 pipeline###
                    if bqsr_dir and vqsr_dir:  ## BQSR and VQSR
                        for sample in path_output:
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample, 'T','T',bqsr_dir))
                        for index, per_cmd in enumerate(cmd_gatk4_call):
                            fw = open(samples[index] + 'ngs_gatk_bqsr_call.sh', 'w').write(per_cmd)
                        files_g_vcf = [i+'.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gcvf.merge(files_g_vcf,'gvcf',out_path+'all_samples',
                                                            [4,ref,genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir])
                        cmd_vqsr = variation_qc.vqsr(ref,out_path+'all_samples.combined.vcf',vqsr_dir,out_path+"all_samples")
                        fw = open('ngs_gatk_bqsr_merge_vqsr.sh','w').write(cmd_merge+'\n'+cmd_vqsr)
                    elif bqsr_dir and not vqsr_dir:  ## BQSR and Hard filtering
                        for sample in path_output:
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample, 'T', 'T', bqsr_dir))
                        for index, per_cmd in enumerate(cmd_gatk4_call):
                            fw = open(samples[index] + 'ngs_gatk_bqsr_call.sh', 'w').write(per_cmd)
                        files_g_vcf = [i + '.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gcvf.merge(files_g_vcf, 'gvcf', out_path + 'all_samples',
                                                    [4, ref, genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir])
                        cmd_hard = variation_qc.hard_filter(ref, out_path + 'all_samples.combined.vcf', out_path + "all_samples")
                        fw = open('ngs_gatk_bqsr_merge_hard_filter.sh', 'w').write(cmd_merge + '\n' + cmd_hard)
                    elif not bqsr_dir and not vqsr_dir:  ## Hard filtering only
                        for sample in path_output:
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample, 'F', 'T', bqsr_dir))
                        for index, per_cmd in enumerate(cmd_gatk4_call):
                            fw = open(samples[index] + 'ngs_gatk_call.sh', 'w').write(per_cmd)
                        files_g_vcf = [i + '.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gcvf.merge(files_g_vcf, 'gvcf', out_path + 'all_samples',
                                                         [4, ref, genomicsdb, chr_list, 1, 'map_file.list', 1, 26,
                                                          tmp_dir])
                        cmd_hard = variation_qc.hard_filter(ref, out_path + 'all_samples.combined.vcf',
                                                            out_path + "all_samples")
                        fw = open('ngs_gatk_merge_hard_filter.sh', 'w').write(cmd_merge + '\n' + cmd_hard)
                    ###combine pipeline###
                    cmd_combine = ngs_vars.samtool_gatk_combine(path_output[index])
                    fw = open('ngs_samtools_gatk_bqsr_vqsr_combine.sh', 'w').write(cmd_combine)

            elif mode == 'SV':
                if callpipe == 'breakdancer':
                    sv_cmd = []
                    for index,sample in path_input:
                        sv_cmd.append(ngs_vars.ngs_sv(sample,'',ref,'breakdancer','false'))
                    for index,per_cmd in enumerate(sv_cmd):
                        fw = open(samples[index] + 'ngs_bd_sv.sh', 'w').write(per_cmd)
                elif callpipe == 'crest':
                    sv_cmd = []
                    for index,sample in path_input:
                        sv_cmd.append(ngs_vars.ngs_sv(sample,'',ref,'crest','false'))
                    for index, per_cmd in enumerate(sv_cmd):
                        fw = open(samples[index] + 'ngs_crest_sv.sh', 'w').write(per_cmd)

            elif mode == 'CNV':
                if callpipe == 'cnvnator':
                    cnv_cmd = []
                    for index, sample in path_input:
                        cnv_cmd.append(ngs_vars.ngs_cnv(sample,'',ref,out_path,'cnvnator','human',strategy,'false'))
                    for index, per_cmd in enumerate(cnv_cmd):
                        fw = open(samples[index] + 'ngs_cnvnator.sh','w').write(per_cmd)
                elif callpipe == 'control-freec':
                    for index, sample in path_input:
                        cnv_cmd.append(ngs_vars.ngs_cnv(sample,'',ref,out_path,'control-freec','human',strategy,'false'))
                    for index, per_cmd in enumerate(cnv_cmd):
                        fw = open(samples[index] + 'ngs_control_freec.sh','w').write(per_cmd)

            elif mode == 'SNP_Indel_Somatic':
                ### gatk4 pipeline ###
                cmd_create_pon = []
                cmd_somatic = ''
                if pon and not normal_samples_for_pon:
                    pon_dir = inputs_dir+pon
                    if germline:
                        cmd_somatic = somatic_detection.mutect2(input_path, out_path, ref, sample,
                                                                                control,
                                                                                pon_dir, germline, af, 3,
                                                                                normal_samples_for_pon_dir)[1]
                    elif not germline:
                        cmd_somatic = somatic_detection.mutect2(input_path, out_path, ref, sample,
                                                                                control,
                                                                                pon_dir, '', '', 3,
                                                                                normal_samples_for_pon_dir)[1]
                    fw = open('somatic_snp_indel.sh','w').write(cmd_somatic)

                elif not pon and normal_samples_for_pon:
                    normal_samples_for_pon_dir = normal_samples_for_pon.strip().split(',')
                    if germline:
                        cmd_create_pon,cmd_somatic = somatic_detection.mutect2(input_path, out_path, ref, sample,
                                                                                control,
                                                                                pon_dir, germline, af, 3,
                                                                                normal_samples_for_pon_dir)
                    elif not germline:
                        cmd_create_pon,cmd_somatic = somatic_detection.mutect2(input_path, out_path, ref, sample,
                                                                                control,
                                                                                pon_dir, '', '', 3,
                                                                                normal_samples_for_pon_dir)
                    for index,sample in normal_samples_for_pon_dir:
                        fw = open(sample+'create_pon.sh','w').write(cmd_create_pon[index])
                    fw = open('somatic_snp_indel.sh','w').write(cmd_somatic)


            elif mode == 'SV_Somatic':
                ### crest pipeline###
                sv_cmd = ''
                sv_cmd = ngs_vars.ngs_sv(sample,control,ref,'crest','true')
                fw = open('somatic_sv.sh','w').write(sv_cmd)

            elif mode == 'CNV_Somatic':
                ### control-freec pipeline ###
                snv_cmd = ''
                snv_cmd = ngs_vars.ngs_cnv(sample,control,ref,'control-freec','human',strategy,'true',out_path)

        elif platform == 'tgs':
            if mode == 'SNP_Indel':
                ### gatk4 pipeline for CCS reads###
                snp_indel_cmd = tgs_vars.tgs_snp_indel(ref,)
            elif mode == 'SV':
                ### sniffles pipeline ###

def write_annotation():
    try:
        input_path = os.path.abspath(outputs_dir) + '/var_results/'
        out_path = os.path.abspath(outputs_dir) + '/annotation_results/'
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
    mutation.add_argument('-sg', '--strategy', type=str, default='WGS',choices=['WGS','WES'],
                          help='WGS or WES')
    mutation.add_argument('-cp','--callpipe',type=str,default='samtools',choices=['samtools','gatk4','samtools+gakt4',
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
    mutation.add_argument('-bv','--build_version',type=str,default='hg19',choices=['hg19','hg38','other'],
                          help='Species and genome version')
    mutation.add_argument('-vc','--variation_calling', type=str,default='single',choices=['single','join'],
                           help='Calling a group of samples together(join calling) or Variant calling with a single \
                           sample only(single sample calling), default is single')  ### https://bcbio.wordpress.com/2014/10/07/joint-calling/
    mutation.add_argument('-s1', '--sample',type=str,
                          help='The target(tumor/disease) sample name')
    mutation.add_argument('-s2', '--control', type=str,
                           help='The control(normal) sample name')
    mutation.add_argument('-pon', '--panel_of_normals', type=str,default='',
                          help='The PON(panel of normal) vcf file,if not provided, use the --normal_samples_for_pon paramters \
                               to create the normal panel')
    mutation.add_argument('-np', '--normal_samples_for_pon', type=str,default='',
                          help='The sample name of normal samples to create the PON,e.g. sample1,sample2,sample3...')
    mutation.add_argument('-gm', '--germline', type=str,default='',
                          help='The population germline resource')
    mutation.add_argument('-af', '--af_of_alleles_not_in_resource', type=float,default=0.001,
                          help='The GATK4 af-of-alleles-not-in-resource parameter in Mutect2,\
                          The default of 0.001 is appropriate for human sample analyses without any population resource.'
                               )

    args = parser.parse_args()
    inputs_dir = args.input
    outputs_dir = args.outputs
    ref = args.reference
    gff3 = args.gff3
    platform = args.seq_platform
    seqtype = args.seq_type
    maptool = args.maptool
    maptool_parameters = args.maptool_parameters
    callpipe = args.callpipe
    mode = args.mode
    bqsr_dir = args.BQSR
    vqsr_dir = args.VQSR
    v_calling = args.variation_calling
    sample = args.sample
    control = args.control
    strategy = args.strategy
    buildver = args.build_version
    pon = args.panel_of_normals
    normal_samples_for_pon = args.normal_samples_for_pon
    germline = args.germline
    af = args.af_of_alleles_not_in_resource

    write_mapping()
    write_call_var()
    write_annotation()



