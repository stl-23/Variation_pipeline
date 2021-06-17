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
    subprocess.check_call(['mkdir', '-p', '{0}'.format(out_path)])
    # lst = os.listdir(input_path)
    #outfile = []
    global samples
    samples = set()

    if platform == "ngs":
        samples = parsering.parse_short_read_dir(input_path, out_path, seqtype)[0]
        samples = list(samples)
        outfile = mapping.Ngs(maptool, input_path, out_path,prefix,maptool_parameters,seqtype).ngs()
        for i,s in enumerate(outfile):
            with open('s1_ngs_'+samples[i]+'.bwa.sh','w') as fw:
                fw.write(s)

    elif platform == 'tgs':
        samples = parsering.parse_long_read_dir(input_path)[0]
        samples = list(samples)
        if maptool == 'Minimap2':
            outfile = mapping.Tgs(maptool,input_path,out_path,ref,maptool_parameters).tgs_minimap2()
            for i,s in enumerate(outfile):
                with open('s1_tgs_'+samples[i]+'.minimap2.sh','w') as fw:
                    fw.write(s)
        elif maptool == 'NGMLR':
            outfile = mapping.Tgs(maptool,input_path,out_path,ref,maptool_parameters).tgs_ngmlr()
            for i,s in enumerate(outfile):
                with open('s1_tgs_'+samples[i]+'.ngmlr.sh','w') as fw:
                    fw.write(s)


def write_call_var():
    input_path = os.path.join(os.path.abspath(outputs_dir), 'mapping_results/')
    out_path = os.path.join(os.path.abspath(outputs_dir), 'var_results/')
    tmp_dir = os.path.join(os.path.abspath(outputs_dir), 'tmp_dir/')
    subprocess.check_call(['mkdir','-p',out_path,tmp_dir])
    # lst = os.listdir(input_path)
    path_input = [input_path + x+'.rmdup.bam' for x in samples]
    path_output = [out_path + x for x in samples]
    map_file = []

    if platform == 'ngs':
        if mode == 'SNP_INDEL':
            for s in samples:
                map_file.append(s + '\t' + s + '.g.vcf')
            with open('map_file.list', 'w') as fw:
                fw.write('\n'.join(map_file) + '\n')
            if callpipe == 'samtools':
                cmd_call = []
                cmd_merge = ''
                if v_calling == 'single': ## single sample calling
                    for index,sample in enumerate(path_output):
                        cmd_call.append(ngs_vars.snp_indel_samtools(ref, path_input[index], sample, v_calling,
                                                                    bcftools_filter))
                    for index, per_cmd in enumerate(cmd_call):
                        with open('s2_ngs_'+samples[index]+'_samtools_call.sh','w') as fw:
                            fw.write(per_cmd)
                elif v_calling == 'join': ## join calling
                    out_name = out_path+'all_samples'
                    cmd_call.append(ngs_vars.snp_indel_samtools(ref, path_input, out_name, v_calling,
                                                                bcftools_filter
                                                                ))
                    with open('s2_ngs_all_samples_samtools_call.sh','w') as fw:
                        fw.write(cmd_call[0])
            elif callpipe == 'gatk4':
                if bqsr_dir and vqsr_dir: ## BQSR and VQSR
                    cmd_call = []
                    cmd_merge = ''
                    if v_calling == 'single': ##  single sample calling
                        print("Warning: too few samples, VQSR may be failed!")
                        cmd_vqsr = []
                        for index,sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample,'F',bqsr_dir))
                            cmd_vqsr.append(variation_qc.vqsr(ref,path_output[index]+'.vcf',vqsr_dir,path_output[index])) ### VQSR ###
                        for i, s in enumerate(samples):
                            with open('s2_ngs_'+s+'_gatk_bqsr_call_vqsr.sh','w') as fw:
                                fw.write(cmd_call[i]+'\n'+cmd_vqsr[i]+'\n')
                    elif v_calling == 'join':  ## join calling
                        cmd_vqsr = ''
                        for index, sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample,'T',bqsr_dir))
                        for index, per_cmd in enumerate(cmd_call):
                            with open('s2.1_ngs_'+samples[index] + '_gatk_bqsr_call.sh', 'w') as fw:
                                fw.write(per_cmd)
                        files_g_vcf = [i+'.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf,'gvcf',out_path+'all_samples',
                                                            maxm,ref,genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir)
                        cmd_vqsr = variation_qc.vqsr(ref,out_path+'all_samples.combined.vcf',vqsr_dir,out_path+"all_samples")
                        with open('s2.2_ngs_gatk_bqsr_merge_vqsr.sh','w') as fw:
                            fw.write(cmd_merge+'\n'+cmd_vqsr)

                elif bqsr_dir and not vqsr_dir: ## BQSR and Hard filtering
                    cmd_call = []
                    cmd_merge = ''
                    if v_calling == 'single':  ##  single sample calling
                        cmd_hard = []
                        for index, sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'F', bqsr_dir))
                            cmd_hard.append(variation_qc.hard_filter(ref, path_output[index]+'.vcf', path_output[index]))
                        for i, s in enumerate(samples):
                            with open('s2_ngs_'+s + '_gatk_bqsr_call_hard_filter.sh', 'w') as fw:
                                fw.write(cmd_call[i] + '\n' + cmd_hard[i] + '\n')
                    elif v_calling == 'join':  ## join calling
                        cmd_hard = ''
                        for index, sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'T', bqsr_dir))
                        for index, per_cmd in enumerate(cmd_call):
                            with open('s2.1_ngs_'+samples[index] + '_gatk_bqsr_call.sh', 'w') as fw:
                                fw.write(per_cmd)
                        files_g_vcf = [i + '.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf, 'gvcf', out_path + 'all_samples',
                                                    maxm, ref, genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir)
                        cmd_hard = variation_qc.hard_filter(ref, out_path + 'all_samples.combined.vcf', out_path + "all_samples")
                        with open('s2.2_ngs_gatk_bqsr_merge_hard_filter.sh', 'w') as fw:
                            fw.write(cmd_merge + '\n' + cmd_hard)

                elif not bqsr_dir and not vqsr_dir: ## Hard filtering only
                    cmd_call = []
                    cmd_merge = ''
                    if v_calling == 'single':  ##  single sample calling
                        cmd_hard = []
                        for index, sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample,'F', bqsr_dir))
                            cmd_hard.append(variation_qc.hard_filter(ref, path_output[index] + '.vcf', path_output[index]))
                        for i, s in enumerate(samples):
                            with open('s2_ngs_'+s + '_gatk_call_hard_filter.sh', 'w') as fw:
                                fw.write(cmd_call[i] + '\n' + cmd_hard[i] + '\n')
                    elif v_calling == 'join':  ## join calling
                        cmd_hard = ''
                        for index, sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'T', bqsr_dir))
                        for index, per_cmd in enumerate(cmd_call):
                            with open('s2.1_ngs_'+samples[index] + '_gatk_call.sh', 'w') as fw:
                                fw.write(per_cmd)
                        files_g_vcf = [i + '.g.vcf' for i in path_output]

                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf, 'gvcf', out_path + 'all_samples',
                                                         maxm, ref, genomicsdb, chr_list, 1, 'map_file.list', 1, 26,tmp_dir)
                        cmd_hard = variation_qc.hard_filter(ref, out_path + 'all_samples.combined.vcf',
                                                            out_path + "all_samples")
                        with open('s2.2_ngs_gatk_merge_hard_filter.sh', 'w') as fw:
                            fw.write(cmd_merge + '\n' + cmd_hard)
                elif not bqsr_dir and vqsr_dir:
                    cmd_call = []
                    cmd_merge = ''
                    if v_calling == 'single': ##  single sample calling
                        print("Warning: too few samples, VQSR may be failed!")
                        cmd_vqsr = []
                        for index,sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample,'F',bqsr_dir))
                            cmd_vqsr.append(variation_qc.vqsr(ref,path_output[index]+'.vcf',vqsr_dir,path_output[index])) ### VQSR ###
                        for i, s in enumerate(samples):
                            with open('s2_ngs_'+s+'_gatk_call_vqsr.sh','w') as fw:
                                fw.write(cmd_call[i]+'\n'+cmd_vqsr[i]+'\n')
                    elif v_calling == 'join':  ## join calling
                        cmd_vqsr = ''
                        for index, sample in enumerate(path_output):
                            cmd_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample,'T',bqsr_dir))
                        for index, per_cmd in enumerate(cmd_call):
                            with open('s2.1_ngs_'+samples[index] + '_gatk_call.sh', 'w') as fw:
                                fw.write(per_cmd)
                        files_g_vcf = [i+'.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf,'gvcf',out_path+'all_samples',
                                                            maxm,ref,genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir)
                        cmd_vqsr = variation_qc.vqsr(ref,out_path+'all_samples.combined.vcf',vqsr_dir,out_path+"all_samples")
                        with open('s2.2_ngs_gatk_merge_vqsr.sh','w') as fw:
                            fw.write(cmd_merge+'\n'+cmd_vqsr)



            elif callpipe == 'samtools+gatk4':
                if v_calling == 'single':  ##  single sample calling
                    cmd_samtools_call = []
                    cmd_gatk4_call = []
                    cmd_combine = []
                    cmd_vqsr = []
                    cmd_hard = []
                    ###samtools pipeline###
                    for index,sample in enumerate(path_output):
                        cmd_samtools_call.append(ngs_vars.snp_indel_samtools(ref, path_input[index], sample,v_calling,
                                                                             bcftools_filter
                                                                             ))
                    for index, per_cmd in enumerate(cmd_samtools_call):
                        with open('s2.1_ngs_'+samples[index]+'_samtools_call.sh','w') as fw:
                            fw.write(per_cmd)
                    ###gatk4 pipeline###
                    if bqsr_dir and vqsr_dir:  ## BQSR and VQSR
                        for index, sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample,'F',bqsr_dir))
                            cmd_vqsr.append(variation_qc.vqsr(ref,path_output[index]+'.vcf',vqsr_dir,path_output[index]))
                        for i, s in enumerate(samples):
                            with open('s2.1_ngs_'+s + '_gatk_bqsr_call_vqsr.sh', 'w') as fw:
                                fw.write(cmd_gatk4_call[i] + '\n' + cmd_vqsr[i] + '\n')
                    elif bqsr_dir and not vqsr_dir: ## BQSR and Hard filtering
                        for index, sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'F', bqsr_dir))
                            cmd_hard.append(variation_qc.hard_filter(ref, path_output[index] + '.vcf', path_output[index]))
                        for i, s in enumerate(samples):
                            with open('s2.1_ngs_'+s + '_gatk_call_hard_filter.sh', 'w') as fw:
                                fw.write(cmd_gatk4_call[i] + '\n' + cmd_hard[i] + '\n')
                    elif not bqsr_dir and not vqsr_dir: ## Hard filtering only
                        for index, sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'F', bqsr_dir))
                            cmd_hard.append(variation_qc.hard_filter(ref, path_output[index] + '.vcf', path_output[index]))
                        for i, s in enumerate(samples):
                            with open('s2.1_ngs_'+s + '_gatk_call_hard_filter.sh', 'w') as fw:
                                fw.write(cmd_gatk4_call[i] + '\n' + cmd_hard[i] + '\n')
                    elif not bqsr_dir and vqsr_dir:
                        for index,sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample,'F',bqsr_dir))
                            cmd_vqsr.append(variation_qc.vqsr(ref,path_output[index]+'.vcf',vqsr_dir,path_output[index])) ### VQSR ###
                        for i, s in enumerate(samples):
                            with open('s2_ngs_'+s+'_gatk_call_vqsr.sh','w') as fw:
                                fw.write(cmd_call[i]+'\n'+cmd_vqsr[i]+'\n')
                    ###combine pipeline###
                    for index, sample in enumerate(path_output):
                        cmd_combine.append(ngs_vars.samtool_gatk_combine(path_output[index]))

                    for i, s in enumerate(samples):
                        with open('s2.2_ngs_'+s + '_combine.sh', 'w') as fw:
                            fw.write(cmd_combine[i] + '\n')

                elif v_calling == 'join':  ## join calling
                    cmd_samtools_call = []
                    cmd_gatk4_call = []
                    cmd_merge = ''
                    cmd_combine = ''
                    cmd_vqsr = ''
                    cmd_hard = ''
                    ###samtools pipeline###
                    out_name = out_path + 'all_samples'
                    cmd_samtools_call.append(ngs_vars.snp_indel_samtools(ref, path_input, out_name, v_calling,
                                                                         bcftools_filter))
                    with open('s2.1_ngs_all_samples_samtools_call.sh', 'w') as fw:
                        fw.write(cmd_samtools_call[0])
                    ###gatk4 pipeline###
                    if bqsr_dir and vqsr_dir:  ## BQSR and VQSR
                        for index, sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'T',bqsr_dir))
                        for index, per_cmd in enumerate(cmd_gatk4_call):
                            with open('s2.1_ngs_'+samples[index] + '_gatk_bqsr_call.sh', 'w') as fw:
                                fw.write(per_cmd)
                        files_g_vcf = [i+'.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf,'gvcf',out_path+'all_samples',
                                                            maxm,ref,genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir)
                        cmd_vqsr = variation_qc.vqsr(ref,out_path+'all_samples.combined.vcf',vqsr_dir,out_path+"all_samples")
                        with open('s2.2_ngs_gatk_bqsr_merge_vqsr.sh','w') as fw:
                            fw.write(cmd_merge+'\n'+cmd_vqsr)
                    elif bqsr_dir and not vqsr_dir:  ## BQSR and Hard filtering
                        for index, sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input, sample, 'T', bqsr_dir))
                        for index, per_cmd in enumerate(cmd_gatk4_call):
                            with open('s2.1_ngs_'+samples[index] + '_gatk_bqsr_call.sh', 'w') as fw:
                                fw.write(per_cmd)
                        files_g_vcf = [i + '.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf, 'gvcf', out_path + 'all_samples',
                                                    maxm, ref, genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir)
                        cmd_hard = variation_qc.hard_filter(ref, out_path + 'all_samples.combined.vcf', out_path + "all_samples")
                        with open('s2.2_ngs_gatk_bqsr_merge_hard_filter.sh', 'w') as fw:
                            fw.write(cmd_merge + '\n' + cmd_hard)
                    elif not bqsr_dir and not vqsr_dir:  ## Hard filtering only
                        for index, sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample, 'T', bqsr_dir))
                        for index, per_cmd in enumerate(cmd_gatk4_call):
                            with open('s2.1_ngs_'+samples[index] + '_gatk_call.sh', 'w') as fw:
                                fw.write(per_cmd)
                        files_g_vcf = [i + '.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf, 'gvcf', out_path + 'all_samples',
                                                         maxm, ref, genomicsdb, chr_list, 1, 'map_file.list', 1, 26,
                                                          tmp_dir)
                        cmd_hard = variation_qc.hard_filter(ref, out_path + 'all_samples.combined.vcf',
                                                            out_path + "all_samples")
                        with open('s2.2_ngs_gatk_merge_hard_filter.sh', 'w') as fw:
                            fw.write(cmd_merge + '\n' + cmd_hard)
                    elif not bqsr_dir and vqsr_dir:
                        for index, sample in enumerate(path_output):
                            cmd_gatk4_call.append(ngs_vars.snp_indel_gatk(ref, path_input[index], sample,'T',bqsr_dir))
                        for index, per_cmd in enumerate(cmd_gatk4_call):
                            with open('s2.1_ngs_'+samples[index] + '_gatk_call.sh', 'w') as fw:
                                fw.write(per_cmd)
                        files_g_vcf = [i+'.g.vcf' for i in path_output]
                        cmd_merge = merge_vcf_gvcf.merge(files_g_vcf,'gvcf',out_path+'all_samples',
                                                            maxm,ref,genomicsdb, chr_list, 1, 'map_file.list', 1, 26, tmp_dir)
                        cmd_vqsr = variation_qc.vqsr(ref,out_path+'all_samples.combined.vcf',vqsr_dir,out_path+"all_samples")
                        with open('s2.2_ngs_gatk_merge_vqsr.sh','w') as fw:
                            fw.write(cmd_merge+'\n'+cmd_vqsr)
                    ###combine pipeline###
                    cmd_combine = ngs_vars.samtool_gatk_combine(out_path+'all_samples')
                    with open('s2.3_ngs_combine.sh', 'w') as fw:
                        fw.write(cmd_combine)


        elif mode == 'SV':
            if callpipe == 'breakdancer':
                sv_cmd = []
                for index,sample in enumerate(path_input):
                    sv_cmd.append(ngs_vars.ngs_sv(sample,'',ref,'breakdancer','false'))
                for index,per_cmd in enumerate(sv_cmd):
                    with open('s2_ngs_'+samples[index] + '_breakdancer_sv.sh', 'w') as fw:
                        fw.write(per_cmd)
            elif callpipe == 'crest':
                sv_cmd = []
                for index,sample in enumerate(path_input):
                    sv_cmd.append(ngs_vars.ngs_sv(sample,'',ref,'crest','false'))
                for index, per_cmd in enumerate(sv_cmd):
                    with open('s2_ngs_'+samples[index] + '_crest_sv.sh', 'w') as fw:
                        fw.write(per_cmd)

        elif mode == 'CNV':
            if callpipe == 'cnvnator':
                cnv_cmd = []
                for index, sample in enumerate(path_input):
                    cnv_cmd.append(ngs_vars.ngs_cnv(sample,'',ref,out_path,'cnvnator','human',strategy,'false'))
                for index, per_cmd in enumerate(cnv_cmd):
                    with open('s2_ngs_'+samples[index] + '_cnvnator.sh','w') as fw:
                        fw.write(per_cmd)
            elif callpipe == 'control-freec':
                cnv_cmd = []
                for sample in path_input:
                    pre = sample.rstrip('/').split('/')[-1]
                    make_conf = ngs_vars.ngs_cnv(sample,'',ref,out_path,'control-freec','human',strategy,'false')[0]
                    with open(pre + '_config_wgs_no_control.list', 'w') as fw:
                        fw.write(make_conf)
                    cnv_cmd.append(ngs_vars.ngs_cnv(sample,'',ref,out_path,'control-freec','human',strategy,'false')[1])
                for index, per_cmd in enumerate(cnv_cmd):
                    with open('s2_ngs_'+samples[index] + '_control_freec.sh','w') as fw:
                        fw.write(per_cmd)

        elif mode == 'SNP_INDEL_Somatic':
            ### gatk4 pipeline ###

            cmd_create_pon = []
            cmd_somatic = ''
            if pon and not normal_samples_for_pon:
                cmd_somatic = somatic_detection.mutect2(input_path, out_path, ref, tar,con,interval_list,pon,
                                                                germline, af, maxm,None)[1]
                with open('s2_somatic_snp.sh','w') as fw:
                    fw.write(cmd_somatic)

            elif not pon and normal_samples_for_pon:

                normal_samples_for_pon_list = normal_samples_for_pon.strip().split(',')
                cmd_create_pon,cmd_somatic = somatic_detection.mutect2(input_path, out_path, ref, tar,con,interval_list,
                                                                        None, germline, af, maxm, normal_samples_for_pon)
                for index,sample in enumerate(normal_samples_for_pon_list):
                    with open('s2.1_'+sample+'_create_pon.sh','w') as fw:
                        fw.write(cmd_create_pon[index])
                with open('s2.2_somatic_snp.sh','w') as fw:
                    fw.write(cmd_somatic)
            elif not pon and not normal_samples_for_pon:
                cmd_somatic = somatic_detection.mutect2(input_path, out_path, ref, tar, con, interval_list, None,
                                                            germline, af, maxm, None)[1]
                with open('s2_somatic_snp.sh', 'w') as fw:
                    fw.write(cmd_somatic)
            elif pon and normal_samples_for_pon:
                print('Warning: parameters -pon and -np were used at the same time, prefer use the PON file first')
                cmd_somatic = somatic_detection.mutect2(input_path, out_path, ref, tar,con,interval_list,pon,
                                                            germline, af, maxm, None)[1]
                with open('s2_somatic_snp.sh','w') as fw:
                    fw.write(cmd_somatic)

        elif mode == 'SV_Somatic':
            ### crest pipeline###
            target = input_path + '/' + tar
            control = input_path + '/' + con
            sv_cmd = ngs_vars.ngs_sv(target,control,ref,'crest','true')
            with open('s2_somatic_sv.sh','w') as fw:
                fw.write(sv_cmd)

        elif mode == 'CNV_Somatic':
            ### control-freec pipeline ###
            target = input_path + '/' + tar
            control = input_path + '/' + con
            config, cnv_cmd = ngs_vars.ngs_cnv(target,control,ref,out_path,'control-freec','human',strategy,'true')
            with open('config.list','w') as fc:
                fc.write(config)
            with open('s2.2_somatic_cnv.sh','w') as fw:
                fw.write(cnv_cmd)

    elif platform == 'tgs':
        if mode == 'SNP_INDEL':
            ### gatk4 pipeline for CCS reads###
            snp_indel_cmds = []
            for index,sample in enumerate(path_input):
                snp_indel_cmds.append(tgs_vars.tgs_snp_indel(ref,sample,path_output[index]))
            for i,s in enumerate(samples):
                with open('s2_tgs_'+s+'_snp_indel_gatk4.sh', 'w') as fw:
                    fw.write(snp_indel_cmds[i])

        elif mode == 'SV':
            ### sniffles pipeline ###
            sv_cmds = []
            for index,sample in enumerate(path_input):
                sv_cmds.append(tgs_vars.tgs_sv(sample,'sniffles',sniffles_p))
            for i,s in enumerate(samples):
                with open('s2_tgs_'+s+'_sv_sniffles.sh','w') as fw:
                    fw.write(sv_cmds[i])

def write_annotation():
    input_path = os.path.join(os.path.abspath(outputs_dir), 'var_results/')
    out_path = os.path.join(os.path.abspath(outputs_dir),'annotation_results/')
    subprocess.run(['mkdir', '-p', out_path])
    #anno_cmd = {}
    if v_calling == 'single':
        for index,sample in enumerate(samples):
            with open('s3_' + sample + '_annotation.sh', 'w') as fw:
                if mode == 'SNP_INDEL' or mode == 'SNP_INDEL_Somatic':
                    snp_cmd = annotation.annotation('annovar', ref, input_path + sample + '.final.snp.gz', gff3,
                                                    out_path + sample + '.final.snp', buildver)
                    indel_cmd = annotation.annotation('annovar', ref, input_path + sample + '.final.indel.gz', gff3,
                                                      out_path + sample + '.final.indel', buildver)
                    fw.write(snp_cmd + '\n' + indel_cmd + '\n')

                elif mode == 'SV' or mode == 'SV_Somatic':
                    fw.write(annotation.annotation('annovar',ref,input_path+sample+'.final.sv.vcf.gz',gff3,
                                                   out_path+sample+'.sv',buildver)+'\n')
        #elif os.path.exists(input_path+sample+'.cnv.vcf'):
                elif mode == 'CNV' or mode == 'CNV_Somatic':
                    fw.write(annotation.annotation('annovar',ref,input_path+sample+'.final.cnv.vcf.gz',gff3,
                                                   out_path+sample+'.cnv',buildver)+'\n')
    elif v_calling == 'join':
        with open('s3_annotation.sh', 'w') as fw:
            if mode == 'SNP_INDEL' or mode == 'SNP_INDEL_Somatic':
                snp_cmd = annotation.annotation('annovar', ref, input_path + 'all_samples.final.snp.gz', gff3,
                                                out_path + 'all_samples.final.snp', buildver)
                indel_cmd = annotation.annotation('annovar', ref, input_path + 'all_samples.final.indel.gz', gff3,
                                                  out_path + 'all_samples.final.indel', buildver)
                fw.write(snp_cmd + '\n' + indel_cmd + '\n')


if __name__ == '__main__':
    examplelog = """EXAMPLES:
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -bv hg38 -sp ngs -mt BWA -cp samtools 
    -sg WGS -mode SNP_INDEL
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -bv hg19 -sp ngs -mt BWA -cp gatk4 
    -bqsr -vqsr -mode SNP_INDEL
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -cp gatk4 -tar Tu4 -con Nm_35 -mode SNP_INDEL
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp ngs -mt BWA -cp control-freec -mode CNV -sg WES 
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -bv hg38 -sp ngs -mt BWA -cp gatk4 
    -mode SNP_INDEL_Somatic -tar SRR1553607.1 -con SRR1553608
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -bv hg38 -sp ngs -mt BWA -cp gatk4 
    -mode SNP_INDEL_Somatic -tar SRR1553607.1 -con SRR1553608 -np sample1,sample2,samples3
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -bv hg38 -sp ngs -mt BWA -cp gatk4 
    -mode SNP_INDEL_Somatic -tar SRR1553607.1 -con SRR1553608 -pon /root/my_data/genomicsdb/g38/hg38_pon.vcf
    -gm /root/my_data/germlines/hg19/af-only-gnomad.vcf.gz -af 0.001 
    python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp tgs -mt Minimap2 -cp gatk4 -mode SNP_INDEL
     python run_variation.py -i /root/my_data/cleandata/ -o /root/my_data/results/ -r /root/my_data/ref/hg19.fa -g 
    /root/my_data/ref/hg19.gff3 -sp tgs -mt NGMLR -cp sniffles -mode SV
    """
    parser = argparse.ArgumentParser(description='Variation pipeline v1.0',
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
    general.add_argument('-mx','--max_memory',type=int,default=24,
                         help='The maximum JAVA memory in the pipeline')
    mapref = parser.add_argument_group(title='Mapping options')
    mapref.add_argument('-sp', '--seq_platform',type=str, default='ngs',choices=['ngs','tgs'],
                        help='Reads sequencing platform, ngs for illumina short reads; tgs for Pacbio or ONT long reads')
    mapref.add_argument('-st', '--seq_type', type=str, default='PE', choices=['PE', 'SE'],
                        help='Sequencing type, paired-end/mate-pair, or single end')
    mapref.add_argument('-mt', '--maptool', type=str, default='BWA', choices=['BWA','Minimap2','NGMLR'],
                        help='Choose an alignment tool,illumina:BWA;Pacbio/ONT:Minimap2/NGMLR')
    mapref.add_argument('-mp', '--maptool_parameters',type=str,default='',
                        help="Set parameters for alignment tools")
    mutation = parser.add_argument_group(title='Detection variation (SNP/INDEL/SV/CNV) options')
    mutation.add_argument('-sg', '--strategy', type=str, default='WGS',choices=['WGS','WES'],
                          help='WGS or WES')
    mutation.add_argument('-cp','--callpipe',type=str,default='samtools',choices=['samtools','gatk4','samtools+gatk4',
                                                                                  'breakdancer','crest',
                                                                                  'cnvnator','control-freec','sniffles'],
                          help='Choose a detection pipeline for SNP/INDEL/SV/CNV')
    mutation.add_argument('-mode', '--mode', type=str, default='SNP_INDEL',choices=['SNP_INDEL','SNP_INDEL_Somatic',
                                                                                    'CNV','SV',
                                                                                    'SV_Somatic','CNV_Somatic'],
                          help='Mutation types')
    mutation.add_argument('--bcftools_filter',type=str,default='',
                          help='Hard filter paramters for SNP/Indel detection pipeline by bcftools. See '
                               'http://samtools.github.io/bcftools/howtos/variant-calling.html')
    mutation.add_argument('-bqsr','--BQSR', nargs='?', const=True, action="store",
                          help='The directory that only contain known sites for BQSR parameter(GATK4 only)')
    mutation.add_argument('-vqsr','--VQSR', nargs='?', const=True, action="store",
                          help='Choose a method for mutations quality control(GATK4 only), if set this parameters, \
                               please provide a directory that only contain resource files; if not, hard_filtering will \
                               be used')
    mutation.add_argument('-vc','--variation_calling', type=str,default='single',choices=['single','join'],
                           help='Calling a group of samples together(join calling) or Variant calling with a single \
                           sample only(single sample calling), default is single')  ### https://bcbio.wordpress.com/2014/10/07/joint-calling/
    somatic = parser.add_argument_group(title='Somatic mutation options')
    somatic.add_argument('-tar', '--target',type=str,
                          help='The target(case/tumor/disease) sample name')
    somatic.add_argument('-con', '--control', type=str,
                           help='The control(normal) sample name')
    somatic.add_argument('-pon', '--panel_of_normals', type=str,default='',
                          help='The PON(panel of normal) vcf file,if not provided, use the --normal_samples_for_pon paramters \
                               to create the normal panel')
    somatic.add_argument('-np', '--normal_samples_for_pon', type=str,default='',
                          help='The sample name of normal samples to create the PON,e.g. sample1,sample2,sample3...')
    somatic.add_argument('-gm', '--germline', type=str,default='',
                          help='The population germline resource')
    somatic.add_argument('-af', '--af_of_alleles_not_in_resource', type=float,default=0.001,
                          help='The GATK4 af-of-alleles-not-in-resource parameter in Mutect2,\
                          The default of 0.001 is appropriate for human sample analyses without any population resource.'
                               )
    somatic.add_argument('--interval',type=str,default='',
                         help='The interval file')
    tgs = parser.add_argument_group(title='TGS(Third generation sequencing) SV options')
    tgs.add_argument('--sniffles_parameters',type=str,default='-s 1 -d 600 --genotype --cluster --ccs_reads',
                          help='The sniffles parameters for SV calling. ')


    args = parser.parse_args()
    inputs_dir = args.input
    outputs_dir = args.output
    ref = args.reference
    gff3 = args.gff3
    platform = args.seq_platform
    seqtype = args.seq_type
    maxm = args.max_memory
    maptool = args.maptool
    maptool_parameters = args.maptool_parameters.strip('"')
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
    sniffles_p = args.sniffles_parameters.strip('"')
    bcftools_filter = args.bcftools_filter.strip('"')
    interval = args.interval

### check arguments ###
    index_shell = os.path.abspath(os.path.dirname(__file__))+'/vartools/index.sh'
    statistics_shell = os.path.abspath(os.path.dirname(__file__))+'/vartools/statistics.sh'
    ## check reference data
    if buildver and not ref and not gff3:
        if buildver == 'hg19':
            genomicsdb = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), 'database/genomicsdb/hg19/')
            ref = os.path.join(genomicsdb,'hg19.fa')
            gff3 = os.path.join(genomicsdb,'hg19.gff3')
            prefix = os.path.splitext(ref)[0]
            chr_list = os.path.join(genomicsdb, 'hg19.chr.list')
            if bqsr_dir:
                bqsr_dir = os.path.join(genomicsdb,'bqsr_resource/')
            if vqsr_dir:
                vqsr_dir = os.path.join(genomicsdb,'vqsr_resource/')
            if not interval:
                interval_list = prefix + '.interval_list'
            if not germline:
                germline = os.path.join(genomicsdb,'pon_and_germline','somatic-hg19_af-only-gnomad.hg19.vcf.gz')
            if not (pon and normal_samples_for_pon):
                pon = os.path.join(genomicsdb,'pon_and_germline','somatic-hg19_1000g_pon.hg19.vcf.gz')
        elif buildver == 'hg38':
            genomicsdb = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), 'database/genomicsdb/hg38/')
            ref = os.path.join(genomicsdb,'hg38.fa')
            prefix = os.path.splitext(ref)[0]
            gff3 = os.path.join(genomicsdb,'hg38.gff3')
            chr_list = os.path.join(os.path.abspath(inputs_dir), 'hg38.chr.list')
            if bqsr_dir:
                bqsr_dir = os.path.join(genomicsdb,'bqsr_resource/')
            if vqsr_dir:
                vqsr_dir = os.path.join(genomicsdb,'vqsr_resource/')
            if not interval:
                interval_list = prefix + '.interval_list'
            if not germline:
                germline = os.path.join(genomicsdb,'pon_and_germline','somatic-hg38_af-only-gnomad.hg38.vcf.gz')
            if not (pon and normal_samples_for_pon):
                pon = os.path.join(genomicsdb,'pon_and_germline','somatic-hg38_1000g_pon.hg38.vcf.gz')
       # else:
       #     genomicsdb = os.path.join(os.path.abspath(inputs_dir), 'genomicsdb/{0}/'.format(buildver))
       #     ref = os.path.join(genomicsdb,'{}.fa'.format(buildver))
       #     gff3 = os.path.join(genomicsdb, '{}.gff'.format(buildver))
       #    chr_list = os.path.join(os.path.abspath(inputs_dir), '{0}_chr.list'.format(buildver))

    elif not buildver and ref and gff3:
        prefix = os.path.splitext(os.path.basename(os.path.abspath(ref)))[0]
        if maptool == 'BWA':
            subprocess.check_call(['sh', index_shell, ref, 'bwa', prefix])
            subprocess.check_call(['sh', index_shell, ref, 'samtools', prefix])
            subprocess.check_call(['sh', statistics_shell, ref, prefix])
            ref = prefix+'.fa'
            gff3 = prefix+'.gff3'
        if mode == 'SNP_INDEL_Somatic' and not interval:
            raise Exception('Error: no interval file provided')

    elif buildver and ref or gff3:
        raise Exception('Do not use -bv and -r/-g together')

    ## check bcftools hard filter paramters
    if bcftools_filter == '':
        if strategy == 'WGS':
            bcftools_filter = """-g3 -G2 -i'%QUAL>=20' -Oz """
        elif strategy == 'WES':
            bcftools_filter = """-sLowQual -g3 -G10 -e'%QUAL<10 || (RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || %MAX(DV)<=3 || %MAX(DV)/%MAX(DP)<=0.3' -Oz"""
    ## check somatic muatation
    if mode == 'SNP_INDEL_Somatic' or mode == 'SV_Somatic' or mode == 'CNV_Somatic':
        if not (tar and con):
            raise Exception('Error: no target and control samples')
    if not germline:
        af = 0

    write_mapping()
    write_call_var()
    write_annotation()



