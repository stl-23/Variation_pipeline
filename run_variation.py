#!/usr/bin/env python
import os
import argparse
import subprocess
import mapping
import ngs_vars
import tgs_vars
import annotation
from vartools import parsering, merge_vcf_gvcf,variation_qc,somatic_detection

def write_mapping():
    input_path = os.path.abspath(inputs_dir) + '/'
    out_path = os.path.join(os.path.abspath(outputs_dir), 'mapping_results/')
    #subprocess.run(['mkdir', '-p', '{0}'.format(out_path)])
    # lst = os.listdir(input_path)
    outfile = []
    global samples
    samples = set()

    if platform == "ngs":
        samples = parsering.parse_short_read_dir(input_path, out_path, seqtype)[0]
        samples = list(samples)
        outfile = mapping.Ngs(maptool, input_path, out_path,ref,maptool_parameters,seqtype).ngs()
        for i,s in enumerate(outfile):
            fw = open('s1_ngs_'+samples[i]+'.bwa.sh','w').write(s)

    elif platform == 'tgs':
        samples = parsering.parse_long_read_dir(input_path)[0]
        samples = list(samples)
        if maptool == 'minimap2':
            outfile = mapping.Tgs(maptool,input_path,out_path,ref,maptool_parameters).tgs_minimap2()
            for i,s in enumerate(outfile):
                fw = open('s1_tgs_'+samples[i]+'.minimap2.sh','w').write(s)
        elif maptool == 'ngml':
            outfile = mapping.Tgs(maptool,input_path,out_path,ref,maptool_parameters).tgs_ngml()
            for i,s in enumerate(outfile):
                fw = open('s1_tgs_'+samples[i]+'.ngml.sh','w').write(s)


def write_call_var():
    input_path = os.path.join(os.path.abspath(outputs_dir), 'mapping_results/')
    out_path = os.path.join(os.path.abspath(outputs_dir), 'var_results/')
    tmp_dir = os.path.join(os.path.abspath(outputs_dir), 'tmp_dir/')
    #subprocess.run(['mkdir','-p',out_path,tmp_dir])
    # lst = os.listdir(input_path)
    path_input = [input_path + x+'.rmdup.bam' for x in samples]
    path_output = [out_path + x for x in samples]
    for s in samples:
        fw = open('map_file.list','a').write(s+'\t'+s+'.g.vcf'+'\n')
    if platform == 'ngs':
        if mode == 'SNP_Indel':
            if callpipe == 'samtools':
                cmd_call = []
                cmd_merge = ''
                if v_calling == 'single': ## single sample calling
                    for index,sample in enumerate(path_output):
                        cmd_call.append(ngs_vars.snp_indel_samtools(ref, path_input[index], sample, v_calling))
                    for index, per_cmd in enumerate(cmd_call):
                        fw = open('s2_ngs_'+samples[index]+'_samtools_call.sh','w').write(per_cmd)
                elif v_calling == 'join': ## join calling
                    out_name = out_path+'all_samples'
                    cmd_call.append(ngs_vars.snp_indel_samtools(ref, path_input, out_name, v_calling))
                    fw = open('s2_ngs_all_samples_samtools_call.sh','w').write(cmd_call[0])
            elif callpipe == 'gatk4':
                if bqsr_dir and vqsr_dir: ## BQSR and VQSR
                    cmd_call = []
                    cmd_merge = ''
                    if v_calling == 'single': ##  single sample calling
                        cmd_vqsr = []
                        for index,sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample,'F',bqsr_dir))
                            cmd_vqsr.append(variation_qc.vqsr(ref,path_output[index]+'.vcf',vqsr_dir,path_output[index]))     ### VQSR ###
                        for i, s in enumerate(samples):
                            fw = open('s2_ngs_'+s+'_gatk_bqsr_call_vqsr.sh','w').write(cmd_call[i]+'\n'+cmd_vqsr[i]+'\n')
                    elif v_calling == 'join':  ## join calling
                        cmd_vqsr = ''
                        for sample in path_output:
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample,'T',bqsr_dir))
                        for index, per_cmd in enumerate(cmd_call):
                            fw = open('s2.1_ngs_'+samples[index] + '_gatk_bqsr_call.sh', 'w').write(per_cmd)
                        files_g_vcf = [i+'.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf,'gvcf',out_path+'all_samples',
                                                            4,ref,genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir)
                        cmd_vqsr = variation_qc.vqsr(ref,out_path+'all_samples.combined.vcf',vqsr_dir,out_path+"all_samples")
                        fw = open('s2.2_ngs_gatk_bqsr_merge_vqsr.sh','w').write(cmd_merge+'\n'+cmd_vqsr)

                elif bqsr_dir and not vqsr_dir: ## BQSR and Hard filtering
                    cmd_call = []
                    cmd_merge = ''
                    if v_calling == 'single':  ##  single sample calling
                        cmd_hard = []
                        for index, sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'F', bqsr_dir))
                            cmd_hard.append(variation_qc.hard_filter(ref, path_output[index]+'.vcf', path_output[index]))
                        for i, s in enumerate(samples):
                            fw = open('s2_ngs_'+s + '_gatk_bqsr_call_hard_filter.sh', 'w').write(cmd_call[i] + '\n' + cmd_hard[i] + '\n')
                    elif v_calling == 'join':  ## join calling
                        cmd_hard = ''
                        for sample in path_output:
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample, 'T', bqsr_dir))
                        for index, per_cmd in enumerate(cmd_call):
                            fw = open('s2.1_ngs_'+samples[index] + '_gatk_bqsr_call.sh', 'w').write(per_cmd)
                        files_g_vcf = [i + '.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf, 'gvcf', out_path + 'all_samples',
                                                    4, ref, genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir)
                        cmd_hard = variation_qc.hard_filter(ref, out_path + 'all_samples.combined.vcf', out_path + "all_samples")
                        fw = open('s2.2_ngs_gatk_bqsr_merge_hard_filter.sh', 'w').write(cmd_merge + '\n' + cmd_hard)

                elif not bqsr_dir and not vqsr_dir: ## Hard filtering only
                    cmd_call = []
                    cmd_merge = ''
                    if v_calling == 'single':  ##  single sample calling
                        cmd_hard = []
                        for index, sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample,'F', bqsr_dir))
                            cmd_hard.append(variation_qc.hard_filter(ref, path_output[index] + '.vcf', path_output[index]))
                        for i, s in enumerate(samples):
                            fw = open('s2_ngs_'+s + '_gatk_call_hard_filter.sh', 'w').write(cmd_call[i] + '\n' + cmd_hard[i] + '\n')
                    elif v_calling == 'join':  ## join calling
                        cmd_hard = ''
                        for sample in path_output:
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample, 'T', bqsr_dir))
                        for index, per_cmd in enumerate(cmd_call):
                            fw = open('s2.1_ngs_'+samples[index] + '_gatk_call.sh', 'w').write(per_cmd)
                        files_g_vcf = [i + '.g.vcf' for i in path_output]

                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf, 'gvcf', out_path + 'all_samples',
                                                         4, ref, genomicsdb, chr_list, 1, 'map_file.list', 1, 26,tmp_dir)
                        cmd_hard = variation_qc.hard_filter(ref, out_path + 'all_samples.combined.vcf',
                                                            out_path + "all_samples")
                        fw = open('s2.2_ngs_gatk_merge_hard_filter.sh', 'w').write(cmd_merge + '\n' + cmd_hard)

            elif callpipe == 'samtools+gatk4':
                if v_calling == 'single':  ##  single sample calling
                    cmd_samtools_call = []
                    cmd_gatk4_call = []
                    cmd_combine = []
                    cmd_vqsr = []
                    cmd_hard = []
                    ###samtools pipeline###
                    for index,sample in enumerate(path_output):
                        cmd_samtools_call.append(ngs_vars.snp_indel_samtools(ref, path_input[index], sample,v_calling))
                    for index, per_cmd in enumerate(cmd_samtools_call):
                        fw = open('s2_ngs_'+samples[index]+'_samtools_call.sh','w').write(per_cmd)
                    ###gatk4 pipeline###
                    if bqsr_dir and vqsr_dir:  ## BQSR and VQSR
                        for index, sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample,'F',bqsr_dir))
                            cmd_vqsr.append(variation_qc.vqsr(ref,path_output[index]+'.vcf',vqsr_dir,path_output[index]))
                        for i, s in enumerate(samples):
                            fw = open('s2_ngs_'+s + '_gatk_bqsr_call_vqsr.sh', 'w').write(cmd_gatk4_call[i] + '\n' + cmd_vqsr[i] + '\n')
                    elif bqsr_dir and not vqsr_dir: ## BQSR and Hard filtering
                        for index, sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'F', bqsr_dir))
                            cmd_hard.append(variation_qc.hard_filter(ref, path_output[index] + '.vcf', path_output[index]))
                        for i, s in enumerate(samples):
                            fw = open('s2_ngs_'+s + '_gatk_call_hard_filter.sh', 'w').write(cmd_gatk4_call[i] + '\n' + cmd_hard[i] + '\n')
                    elif not bqsr_dir and not vqsr_dir: ## Hard filtering only
                        for index, sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'F', bqsr_dir))
                            cmd_hard.append(variation_qc.hard_filter(ref, path_output[index] + '.vcf', path_output[index]))
                        for i, s in enumerate(samples):
                            fw = open('s2_ngs_'+s + '_gatk_call_hard_filter.sh', 'w').write(cmd_gatk4_call[i] + '\n' + cmd_hard[i] + '\n')
                    ###combine pipeline###
                    for index, sample in enumerate(path_output):
                        cmd_combine.append(ngs_vars.samtool_gatk_combine(path_output[index]))
                        for i, s in enumerate(samples):
                            fw = open('s2_ngs_'+s + '_samtools_gatk_bqsr_vqsr_combine.sh', 'w').write(cmd_combine[i] + '\n')

                elif v_calling == 'join':  ## join calling
                    cmd_samtools_call = []
                    cmd_gatk4_call = []
                    cmd_merge = ''
                    cmd_combine = ''
                    cmd_vqsr = ''
                    cmd_hard = ''
                    ###samtools pipeline###
                    out_name = out_path + 'all_samples'
                    cmd_samtools_call.append(ngs_vars.snp_indel_samtools(ref, path_input, out_name, v_calling))
                    fw = open('s2_ngs_all_samples_samtools_call.sh', 'w').write(cmd_samtools_call[0])
                    ###gatk4 pipeline###
                    if bqsr_dir and vqsr_dir:  ## BQSR and VQSR
                        for sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample, 'T',bqsr_dir))
                        for index, per_cmd in enumerate(cmd_gatk4_call):
                            fw = open('s2.1_ngs_'+samples[index] + '_gatk_bqsr_call.sh', 'w').write(per_cmd)
                        files_g_vcf = [i+'.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf,'gvcf',out_path+'all_samples',
                                                            4,ref,genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir)
                        cmd_vqsr = variation_qc.vqsr(ref,out_path+'all_samples.combined.vcf',vqsr_dir,out_path+"all_samples")
                        fw = open('s2.2_ngs_gatk_bqsr_merge_vqsr.sh','w').write(cmd_merge+'\n'+cmd_vqsr)
                    elif bqsr_dir and not vqsr_dir:  ## BQSR and Hard filtering
                        for sample in path_output:
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample, 'T', bqsr_dir))
                        for index, per_cmd in enumerate(cmd_gatk4_call):
                            fw = open('s2.1_ngs_'+samples[index] + '_gatk_bqsr_call.sh', 'w').write(per_cmd)
                        files_g_vcf = [i + '.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf, 'gvcf', out_path + 'all_samples',
                                                    4, ref, genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir)
                        cmd_hard = variation_qc.hard_filter(ref, out_path + 'all_samples.combined.vcf', out_path + "all_samples")
                        fw = open('s2.2_ngs_gatk_bqsr_merge_hard_filter.sh', 'w').write(cmd_merge + '\n' + cmd_hard)
                    elif not bqsr_dir and not vqsr_dir:  ## Hard filtering only
                        for sample in path_output:
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample, 'T', bqsr_dir))
                        for index, per_cmd in enumerate(cmd_gatk4_call):
                            fw = open('s2.1_ngs_'+samples[index] + '_gatk_call.sh', 'w').write(per_cmd)
                        files_g_vcf = [i + '.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf, 'gvcf', out_path + 'all_samples',
                                                         4, ref, genomicsdb, chr_list, 1, 'map_file.list', 1, 26,
                                                          tmp_dir)
                        cmd_hard = variation_qc.hard_filter(ref, out_path + 'all_samples.combined.vcf',
                                                            out_path + "all_samples")
                        fw = open('s2.2_ngs_gatk_merge_hard_filter.sh', 'w').write(cmd_merge + '\n' + cmd_hard)
                    ###combine pipeline###
                        cmd_combine = ngs_vars.samtool_gatk_combine(out_path+'all_samples')
                        fw = open('s2.3_ngs_samtools_gatk_bqsr_vqsr_combine.sh', 'w').write(cmd_combine)

            elif mode == 'SV':
                if callpipe == 'breakdancer':
                    sv_cmd = []
                    for index,sample in enumerate(path_input):
                        sv_cmd.append(ngs_vars.ngs_sv(sample,'',ref,'breakdancer','false'))
                    for index,per_cmd in enumerate(sv_cmd):
                        fw = open('s2_ngs_'+samples[index] + '_bd_sv.sh', 'w').write(per_cmd)
                elif callpipe == 'crest':
                    sv_cmd = []
                    for index,sample in enumerate(path_input):
                        sv_cmd.append(ngs_vars.ngs_sv(sample,'',ref,'crest','false'))
                    for index, per_cmd in enumerate(sv_cmd):
                        fw = open('s2_ngs_'+samples[index] + '_crest_sv.sh', 'w').write(per_cmd)

            elif mode == 'CNV':
                if callpipe == 'cnvnator':
                    cnv_cmd = []
                    for index, sample in enumerate(path_input):
                        cnv_cmd.append(ngs_vars.ngs_cnv(sample,'',ref,out_path,'cnvnator','human',strategy,'false'))
                    for index, per_cmd in enumerate(cnv_cmd):
                        fw = open('s2_ngs_'+samples[index] + '_cnvnator.sh','w').write(per_cmd)
                elif callpipe == 'control-freec':
                    cnv_cmd = []
                    for index, sample in enumerate(path_input):
                        cnv_cmd.append(ngs_vars.ngs_cnv(sample,'',ref,out_path,'control-freec','human',strategy,'false'))
                    for index, per_cmd in enumerate(cnv_cmd):
                        fw = open('s2_ngs_'+samples[index] + '_control_freec.sh','w').write(per_cmd)

            elif mode == 'SNP_Indel_Somatic':
                ### gatk4 pipeline ###
                cmd_create_pon = []
                cmd_somatic = ''
                if pon and not normal_samples_for_pon:
                    pon_dir = inputs_dir+pon
                    if germline:
                        cmd_somatic = somatic_detection.mutect2(input_path, out_path, ref, tar,con,pon_dir,
                                                                germline, af, 3, normal_samples_for_pon)[1]
                    elif not germline:
                        cmd_somatic = somatic_detection.mutect2(input_path, out_path, ref, tar, con, pon_dir, '', '', 3,
                                                                normal_samples_for_pon)[1]
                    fw = open('s2_somatic_snp_indel.sh','w').write(cmd_somatic)

                elif not pon and normal_samples_for_pon:
                    if germline:
                        cmd_create_pon,cmd_somatic = somatic_detection.mutect2(input_path, out_path, ref, tar,con, pon,
                                                                               germline, af, 3, normal_samples_for_pon)
                    elif not germline:
                        cmd_create_pon,cmd_somatic = somatic_detection.mutect2(input_path, out_path, ref, tar, con,
                                                                                pon, '', '', 3, normal_samples_for_pon)
                    for index,sample in normal_samples_for_pon:
                        fw = open('s2_'+sample+'_create_pon.sh','w').write(cmd_create_pon[index])
                    fw = open('s2_somatic_snp_indel.sh','w').write(cmd_somatic)


            elif mode == 'SV_Somatic':
                ### crest pipeline###
                target = input_path + '/' + tar
                control = input_path + '/' + con
                sv_cmd = ngs_vars.ngs_sv(target,control,ref,'crest','true')
                fw = open('s2_somatic_sv.sh','w').write(sv_cmd)

            elif mode == 'CNV_Somatic':
                ### control-freec pipeline ###
                target = input_path + '/' + tar
                control = input_path + '/' + con
                cnv_cmd = ngs_vars.ngs_cnv(target,control,ref,'control-freec','human',strategy,'true',out_path)
                fw = open('s2_somatic_cnv.sh','w').write(cnv_cmd)

        elif platform == 'tgs':
            if mode == 'SNP_Indel':
                ### gatk4 pipeline for CCS reads###
                snp_indel_cmds = []
                for index,sample in enumerate(path_input):
                    snp_indel_cmds.append(tgs_vars.tgs_snp_indel(ref,sample,path_output[index]))
                for i,s in samples:
                    fw = open('s2_tgs_'+s+'_snp_indel_gatk4.sh', 'w').write(snp_indel_cmds[i])

            elif mode == 'SV':
                ### sniffles pipeline ###
                sv_cmds = []
                for index,sample in enumerate(path_input):
                    sv_cmds.append(tgs_vars.tgs_sv(sample,'sniffles',sniffles_p))
                for i,s in samples:
                    fw = open('s2_tgs_'+s+'_sv_sniffles.sh','w').write(sv_cmds[i])

def write_annotation():
    input_path = os.path.abspath(outputs_dir) + '/var_results/'
    out_path = os.path.abspath(outputs_dir) + '/annotation_results/'
    #subprocess.run(['mkdir', '-p', out_path])
    anno_cmd = []
    for index,sample in enumerate(samples):
        if os.path.exists(input_path+sample+'.snp.vcf'):
            anno_cmd.append(annotation.annotation('annovar',ref,input_path+sample+'.snp.vcf',gff3,out_path+sample+'.snp',buildver))
        elif os.path.exists(input_path+sample+'.indel.vcf'):
            anno_cmd.append(annotation.annotation('annovar',ref,input_path+sample+'.indel.vcf',gff3,out_path+sample+'.indel',buildver))
        elif os.path.exists(input_path+sample+'.sv.vcf'):
            anno_cmd.append(annotation.annotation('annovar',ref,input_path+sample+'.sv.vcf',gff3,out_path+sample+'.sv',buildver))
        elif os.path.exists(input_path+sample+'.cnv.vcf'):
            anno_cmd.append(annotation.annotation('annovar',ref,input_path+sample+'.cnv.vcf',gff3,out_path+sample+'.cnv',buildver))

    fw = open('s3_annotation.sh','w').write('\n'.join(anno_cmd))


if __name__ == '__main__':
    examplelog = """EXAMPLES:
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -bv hg38 -sp ngs -mt BWA -cp samtools -mode SNP_Indel
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -bv hg19 -sp ngs -mt BWA -cp gatk4 -bqsr /root/my_data/known_sites/ -vqsr /root/my_data/resources/
    -pr 'N' -tar target_name -con control_name -mode SNP_Indel
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -cp gatk4 -pr 'N' -tar Tu4 -con Nm_35 -mode SNP_Indel
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -cp gatk4 -bqsr /root/my_data/known_sites/ -vqsr /root/my_data/resources/
    -pr 'Y' -tar target_name -con control_name -bv 'hg19' -vc 'join' -np sample1,sample2,sample3 -gm 
    /root/my_data/germlines/hg19/af-only-gnomad.vcf.gz -af 0.001 
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -cp control-freec -pr 'N' -mode CNV -ot WGS 
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp tgs -mt minimap2 -cp sniffles -mode SV 
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp tgs -mt ngml -cp gatk4 -mode SNP_INDEL 
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
    general.add_argument('-g', '--gff3', type=str,default='',
                         help="The gff3 file of reference for variation sites annotation")
    general.add_argument('-bv','--build_version',type=str,default='',choices=['hg19','hg38'],
                          help='Human genome build version,if used,do not set -r and -g')
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
    mutation.add_argument('-mode', '--mode', type=str, default='SNP_Indel',choices=['SNP_Indel','SNP_Indel_Somatic',
                                                                                    'CNV','SV',
                                                                                    'SV_Somatic','CNV_Somatic'],
                          help='Mutation types')
    mutation.add_argument('-bqsr','--BQSR',type=str,default='',
                          help='The directory that only contain known sites for BQSR parameter(GATK4 only)')
    mutation.add_argument('-vqsr','--VQSR',type=str,default='',
                          help='Choose a method for mutations quality control(GATK4 only), if set this parameters, \
                               please provide a directory that only contain resource files; if not, hard_filtering will \
                               be used')
    mutation.add_argument('-vc','--variation_calling', type=str,default='single',choices=['single','join'],
                           help='Calling a group of samples together(join calling) or Variant calling with a single \
                           sample only(single sample calling), default is single')  ### https://bcbio.wordpress.com/2014/10/07/joint-calling/
    mutation.add_argument('-tar', '--target',type=str,
                          help='The target(tumor/disease) sample name')
    mutation.add_argument('-con', '--control', type=str,
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
    mutation.add_argument('-sniffles_p','--sniffles_parameters',type=str,default='-s 1 -d 600 --genotype --cluster --ccs_reads',
                          help='The sniffles parameters for SV calling. ')

    args = parser.parse_args()
    inputs_dir = args.input
    outputs_dir = args.output
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
    tar = args.target
    con = args.control
    strategy = args.strategy
    buildver = args.build_version
    pon = args.panel_of_normals
    normal_samples_for_pon = args.normal_samples_for_pon
    germline = args.germline
    af = args.af_of_alleles_not_in_resource
    sniffles_p = args.sniffles_parameters
### arguments
    index_shell = os.path.abspath(os.path.dirname(__file__))+'/vartools/index.sh'
    statistics_shell = os.path.abspath(os.path.dirname(__file__))+'/vartools/statistics.sh'
    if buildver and not ref and not gff3:
        if buildver == 'hg19':
            genomicsdb = os.path.join(os.path.abspath(inputs_dir), 'genomicsdb/hg19/')
            ref = os.path.join(genomicsdb,'hg19.fa')
            gff3 = os.path.join(genomicsdb,'hg19.gff')
            chr_list = os.path.join(genomicsdb, 'hg19_chr.list')
            if not bqsr_dir:
                bqsr_dir = os.path.join(genomicsdb,'/bqsr_resource/')
            if not vqsr_dir:
                vqsr_dir = os.path.join(genomicsdb,'/vqsr_resource/')
        elif buildver == 'hg38':
            genomicsdb = os.path.join(os.path.abspath(inputs_dir), 'genomicsdb/hg38/')
            ref = os.path.join(genomicsdb,'hg38.fa')
            gff3 = os.path.join(genomicsdb,'hg38.gff')
            chr_list = os.path.join(os.path.abspath(inputs_dir), 'hg38_chr.list')
            if not bqsr_dir:
                bqsr_dir = os.path.join(genomicsdb,'/bqsr_resource/')
            if not vqsr_dir:
                vqsr_dir = os.path.join(genomicsdb,'/vqsr_resource/')
       # else:
       #     genomicsdb = os.path.join(os.path.abspath(inputs_dir), 'genomicsdb/{0}/'.format(buildver))
       #     ref = os.path.join(genomicsdb,'{}.fa'.format(buildver))
       #     gff3 = os.path.join(genomicsdb, '{}.gff'.format(buildver))
       #    chr_list = os.path.join(os.path.abspath(inputs_dir), '{0}_chr.list'.format(buildver))
    elif not buildver and ref and gff3:
        prefix = os.path.splitext(os.path.basename(os.path.abspath(ref)))[0]
        if maptool == 'BWA':
            subprocess.check_call(['sh', index_shell, ref, 'BWA', prefix])
            subprocess.check_call(['sh', index_shell, ref, 'samtools', prefix])
            subprocess.check_call(['sh', statistics_shell, ref, prefix])
            ref = prefix+'.fa'
            gff3 = prefix+'.gff3'

    elif buildver and ref or gff3:
        print('Do not set -bv and -r/-g together')
        exit(0)

    write_mapping()
    write_call_var()
    write_annotation()



