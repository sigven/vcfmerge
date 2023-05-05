#!/usr/bin/env python

import argparse
import re
import os
import logging
import sys
import gzip

version = '0.1.0'
chromosomes = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']

##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>

def __main__():
   parser = argparse.ArgumentParser(description='Merge somatic calls (SNVs/InDels) from multiple VCF files into a single VCF', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('tumor_sample_id', help='Sample ID for the tumor sample')
   parser.add_argument('control_sample_id', help='Sample ID for the control sample')
   parser.add_argument('output_dir', help='Directory for output files')
   parser.add_argument("--mutect_vcf",dest = "mutect_vcf", help='Bgzipped VCF input file with somatic query variants (SNVs) called with MuTect (version 2.x).')
   parser.add_argument("--strelka_snv_vcf", dest = "strelka_snv_vcf", help='Bgzipped VCF input file with somatic query variants (SNVs) called with Strelka (version 2.x).')
   parser.add_argument("--strelka_indel_vcf", dest = "strelka_indel_vcf", help="Bgzipped VCF input file with somatic query variants (InDels) called with Strelka (version 2.x)")
   parser.add_argument("--compress", action="store_true", help="Compress output VCF with bgzip + tabix")
   parser.add_argument("--force_overwrite", action="store_true", help="Overwrite existing output files")
   args = parser.parse_args()
   
   #mutect_vcf = None
   logger = getlogger('vcf_merge')
   ##check the existence of the output directory
   if not os.path.isdir(args.output_dir):
       logger.error('Output directory ' + str(args.output_dir) + ' does not exist')
       exit(0)

   merge_all_vcf = os.path.join(str(args.output_dir), str(args.tumor_sample_id) + '_' + str(args.control_sample_id) + '_vcfmerge_all.vcf')
   merge_somatic_vcf = os.path.join(str(args.output_dir), str(args.tumor_sample_id) + '_' + str(args.control_sample_id) + '_vcfmerge_somatic.vcf')
   if (os.path.exists(merge_all_vcf) or os.path.exists(str(merge_all_vcf) + '.gz')) and args.force_overwrite is False:
       logger.error('Output file ' + str(merge_all_vcf) + '(.gz) exists - turn on \'--force_overwrite\' to overwrite existing output files')
       exit(0)
   if (os.path.exists(merge_somatic_vcf) or os.path.exists(str(merge_somatic_vcf) + '.gz')) and args.force_overwrite is False:
       logger.error('Output file ' + str(merge_somatic_vcf)  + '(.gz) exists - turn on \'--force_overwrite\' to overwrite existing output files')
       exit(0)
   logger.info('Merging Strelka2 and MuTect2 calls for tumor-normal pair ' + str(args.tumor_sample_id) + '_' + str(args.control_sample_id))

   merge_multiple_vcfs(args.tumor_sample_id, args.control_sample_id, merge_all_vcf, merge_somatic_vcf, logger, args.mutect_vcf, args.strelka_snv_vcf, args.strelka_indel_vcf, args.compress)


def getlogger(logger_name):
   logger = logging.getLogger(logger_name)
   logger.setLevel(logging.DEBUG)

   # create console handler and set level to debug
   ch = logging.StreamHandler(sys.stdout)
   ch.setLevel(logging.DEBUG)

   # add ch to logger
   logger.addHandler(ch)
   
   # create formatter
   formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", "20%y-%m-%d %H:%M:%S")
   
   #add formatter to ch
   ch.setFormatter(formatter)
   
   return logger


def get_sample_column_index(header_line, tumor_sample_id, control_sample_id):
    headers = header_line.rstrip().split('\t')
    i = 0
    tumor_sample_index_column = -1
    control_sample_index_column = -1
    while i < len(headers):
        if headers[i] == tumor_sample_id or str(headers[i]) == 'TUMOR':
            tumor_sample_index_column = i
        if headers[i] == control_sample_id or str(headers[i]) == 'NORMAL':
            control_sample_index_column = i
        i = i + 1
    return {'tumor_sample_index_column':tumor_sample_index_column, 'control_sample_index_column':control_sample_index_column}


def get_algo_prefix(algorithm = 'mutect2'):
    prefix = 'MTCT2'
    if algorithm == 'strelka2':
        prefix = 'STKA2'
    return prefix

def init_genotype(algorithm = 'mutect2'):
    genotype = {}
    algo_prefix = get_algo_prefix(algorithm)
    genotype[str(algo_prefix) + '_FILTER'] = '.'
    genotype[str(algo_prefix) + '_AD_TUMOR'] = '.'
    genotype[str(algo_prefix) + '_AD_CONTROL'] = '.'
    genotype[str(algo_prefix) + '_DP_TUMOR'] = '.'
    genotype[str(algo_prefix) + '_DP_CONTROL'] = '.'
    genotype[str(algo_prefix) + '_VAF_CONTROL'] = '.'
    genotype[str(algo_prefix) + '_VAF_TUMOR'] = '.'
    genotype[str(algo_prefix) + '_CALL_STATISTIC'] = '.'
    genotype[str(algo_prefix) + '_FAIL_REASONS'] = '.'
    return(genotype)


def add_meta_info(meta_info, algorithm = 'mutect2'):

    algorithm_meta = '.'
    if algorithm == 'mutect2':
        algorithm_meta = 'MuTect (2.x)'  
    if algorithm == 'strelka2':
        algorithm_meta = 'Strelka (2.x)'

    algo_prefix = get_algo_prefix(algorithm)
    meta_info[str(algo_prefix) + '_FILTER'] = "##INFO=<ID=" + str(algo_prefix) + "_FILTER" + ",Number=.,Type=String,Description=\"VCF variant FILTER - " + str(algorithm_meta) + "\">"
    meta_info[str(algo_prefix) + '_AD_TUMOR'] = "##INFO=<ID=" + str(algo_prefix) + "_AD_TUMOR" + ",Number=.,Type=String,Description=\"Allelic depths in tumor sample (ref, alt) - " + str(algorithm_meta) + "\">"
    meta_info[str(algo_prefix) + '_AD_CONTROL'] = "##INFO=<ID=" + str(algo_prefix) + "_AD_CONTROL" + ",Number=.,Type=String,Description=\"Allelic depths in control sample (ref, alt) - " + str(algorithm_meta) + "\">"
    meta_info[str(algo_prefix) + '_DP_TUMOR'] = "##INFO=<ID=" + str(algo_prefix) + "_DP_TUMOR" + ",Number=1,Type=Integer,Description=\"Sequencing depth at variant site in tumor sample - " + str(algorithm_meta) + "\">"
    meta_info[str(algo_prefix) + '_DP_CONTROL'] = "##INFO=<ID=" + str(algo_prefix) + "_DP_CONTROL" + ",Number=1,Type=Integer,Description=\"Sequencing depth at variant site in control sample - " + str(algorithm_meta) + "\">"
    meta_info[str(algo_prefix) + '_VAF_TUMOR'] = "##INFO=<ID=" + str(algo_prefix) + "_VAF_CONTROL" + ",Number=1,Type=Float,Description=\"Variant allelic fraction somatic (tumor sample) - " + str(algorithm_meta) + "\">"
    meta_info[str(algo_prefix) + '_VAF_CONTROL'] = "##INFO=<ID=" + str(algo_prefix) + "_VAF_TUMOR" + ",Number=1,Type=Float,Description=\"Variant allelic fraction control sample - " + str(algorithm_meta) + "\">"
    if algorithm == 'mutect2':
        meta_info[str(algo_prefix) + '_CALL_STATISTIC'] = "##INFO=<ID=" + str(algo_prefix) + "_CALL_STATISTIC" + ",Number=.,Type=String,Description=\"MuTect's core statistic: Log of likelihood tumor event is real/likelihood event is sequencing error: TLOD - " + str(algorithm_meta) + "\">"
       
    if algorithm == 'strelka2':
        meta_info[str(algo_prefix) + '_CALL_STATISTIC'] = "##INFO=<ID=" + str(algo_prefix) + "_CALL_STATISTIC" + ",Number=.,Type=String,Description=\"Strelka's quality score: SomaticEVS - " + str(algorithm_meta) + "\">"
    meta_info[str(algo_prefix) + '_FAIL_REASONS'] = "##INFO=<ID=" + str(algo_prefix) + "_FAIL_REASONS" + ",Number=.,Type=String,Description=\"Reasons for not classifying given mutation as somatic - " + str(algorithm_meta) + "\">"
    
    if not 'PASS' in meta_info.keys():
        meta_info['PASS'] = "##FILTER=<ID=PASS,Description=\"All filters passed\">"
    if not 'VARIANT_CALLERS' in meta_info.keys():
        meta_info['VARIANT_CALLERS'] = "##INFO=<ID=VARIANT_CALLERS,Number=.,Type=String,Description=\"Somatic mutation callers that called this variant as somatic (PASS)\">"
    if not 'TDP' in meta_info.keys():
        meta_info['TDP'] = "##INFO=<ID=TDP,Number=1,Type=Integer,Description=\"Sequencing depth at variant site in tumor sample (values from MuTect2 calls have priority over Strelka2)\">"
    if not 'CDP' in meta_info.keys():
        meta_info['CDP'] = "##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Sequencing depth at variant site in control sample (values from MuTect2 calls have priority over Strelka2)\">"
    if not 'TVAF' in meta_info.keys():
        meta_info['TVAF'] = "##INFO=<ID=TVAF,Number=1,Type=Float,Description=\"Variant allelic fraction somatic (tumor sample) (values from MuTect2 calls have priority over Strelka2)\">"
    if not 'CVAF' in meta_info.keys():
        meta_info['CVAF'] = "##INFO=<ID=CVAF,Number=1,Type=Float,Description=\"Variant allelic fraction control sample (values from MuTect2 calls have priority over Strelka2)\">"
    if not 'VT' in meta_info.keys():
        meta_info['VT'] = "##INFO=<ID=VT,Number=.,Type=String,Description=\"Variant type (snv/inDEL/INdel/blocksub\">"
    if not 'SOMATIC' in meta_info.keys():
        meta_info['SOMATIC'] = "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic event (by at least one caller)\">"
    if not 'MNV_SUPPORT_STRELKA' in meta_info.keys():
        meta_info['MNV_SUPPORT_STRELKA'] = "##INFO=<ID=MNV_SUPPORT_STRELKA,Number=0,Type=Flag,Description=\"Multiple Strelka SNVs computationally aggregated into an MNV\">"

def add_filter(algorithm, genotype, filter_val):
    algo_prefix = get_algo_prefix(algorithm)
    genotype[str(algo_prefix) + '_FILTER'] = filter_val

def get_info_data(vcf_info_elements):

    info_columns = {}
    for elem in vcf_info_elements:
        if '=' in elem:
            tag, value = elem.split('=')
            info_columns[tag] = value
        else:
            info_columns[elem] = True
    
    return(info_columns)

def get_gt_data(gt_elements, format_tags):

    gt_data = {}
    if len(gt_elements) == len(format_tags):
        i = 0
        while i < len(gt_elements):
            tag = str(format_tags[i])
            val = str(gt_elements[i])
            gt_data[tag] = val
            i = i + 1
    
    return(gt_data)


def is_float_value(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def populate_genotype_info_tags(genotype_data, ref_allele, alt_allele, info_raw, gt_tumor_raw, gt_control_raw, algorithm, vartype = 'snv'):

    algo_prefix = get_algo_prefix(algorithm)
    nreads_control = 0
    nreads_tumor = 0
    nreads_alt_control = 0
    nreads_alt_tumor = 0
    nreads_ref_control = 0
    nreads_ref_tumor = 0
    control_vaf = 0
    tumor_vaf = 0

    if algorithm == 'mutect2':

        nreads_control = int(gt_control_raw['DP'])
        nreads_tumor = int(gt_tumor_raw['DP'])
        tumor_vaf = "."
        control_vaf = "."

        if is_float_value(gt_tumor_raw['AF']):
            tumor_vaf = float(gt_tumor_raw['AF'])
        
        if is_float_value(gt_control_raw['AF']):
            control_vaf = float(gt_control_raw['AF'])

        if "," in gt_tumor_raw['AD']:
            if str(gt_tumor_raw['AD']).count(',') == 1:
                nreads_ref_tumor = int(str(gt_tumor_raw['AD']).split(',')[0])
                nreads_alt_tumor = int(str(gt_tumor_raw['AD']).split(',')[1])
    
        if "," in gt_control_raw['AD']:
            if str(gt_control_raw['AD']).count(',') == 1:
                nreads_ref_control = int(str(gt_control_raw['AD']).split(',')[0])
                nreads_alt_control = int(str(gt_control_raw['AD']).split(',')[1])
            
        if 'AS_FilterStatus' in info_raw.keys():
            genotype_data[str(algo_prefix) + '_FAIL_REASONS'] = str(info_raw['AS_FilterStatus'])
        if 'TLOD' in info_raw.keys():
            genotype_data[str(algo_prefix) + '_CALL_STATISTIC'] = info_raw['TLOD']

    if algorithm == 'strelka2':
        nreads_control = int(gt_control_raw['DP'])
        nreads_tumor = int(gt_tumor_raw['DP'])

        if vartype == 'snv':
            nreads_control = int(gt_control_raw['DP']) - int(gt_control_raw['FDP'])
            nreads_tumor = int(gt_tumor_raw['DP']) - int(gt_tumor_raw['FDP'])
                 
            if 'SomaticEVS' in info_raw.keys():
                genotype_data[str(algo_prefix) + '_CALL_STATISTIC'] = info_raw['SomaticEVS']
            if ref_allele == 'A':
                nreads_ref_control = int(gt_control_raw['AU'].split(',')[0])
                nreads_ref_tumor = int(gt_tumor_raw['AU'].split(',')[0])
            if ref_allele == 'T':
                nreads_ref_control = int(gt_control_raw['TU'].split(',')[0])
                nreads_ref_tumor = int(gt_tumor_raw['TU'].split(',')[0])
            if ref_allele == 'C':
                nreads_ref_control = int(gt_control_raw['CU'].split(',')[0])
                nreads_ref_tumor = int(gt_tumor_raw['CU'].split(',')[0])
            if ref_allele == 'G':
                nreads_ref_control = int(gt_control_raw['GU'].split(',')[0])
                nreads_ref_tumor = int(gt_tumor_raw['GU'].split(',')[0])
            if alt_allele == 'A':
                nreads_alt_tumor = int(gt_tumor_raw['AU'].split(',')[0])
                nreads_alt_control = int(gt_control_raw['AU'].split(',')[0])
            if alt_allele == 'T':
                nreads_alt_tumor = int(gt_tumor_raw['TU'].split(',')[0])
                nreads_alt_control = int(gt_control_raw['TU'].split(',')[0])
            if alt_allele == 'C':
                nreads_alt_tumor = int(gt_tumor_raw['CU'].split(',')[0])
                nreads_alt_control = int(gt_control_raw['CU'].split(',')[0])
            if alt_allele == 'G':
                nreads_alt_tumor = int(gt_tumor_raw['GU'].split(',')[0])
                nreads_alt_control = int(gt_control_raw['GU'].split(',')[0])

        else:
            nreads_alt_tumor = int(gt_tumor_raw['TIR'].split(',')[0])
            nreads_alt_control = int(gt_control_raw['TIR'].split(',')[0])
            nreads_ref_control = int(gt_control_raw['TAR'].split(',')[0])
            nreads_ref_tumor = int(gt_tumor_raw['TAR'].split(',')[0])

            if 'SomaticEVS' in info_raw.keys():
                genotype_data[str(algo_prefix) + '_CALL_STATISTIC'] = info_raw['SomaticEVS']

    if (nreads_control) > 0 and algorithm == "strelka2":
        control_vaf = float(nreads_alt_control) / nreads_control
    if (nreads_tumor) > 0 and algorithm == "strelka2":
        tumor_vaf = float(nreads_alt_tumor) / nreads_tumor


    genotype_data[str(algo_prefix) + '_VAF_TUMOR'] = "{0:.4f}".format(tumor_vaf)
    genotype_data[str(algo_prefix) + '_VAF_CONTROL'] = "{0:.4f}".format(control_vaf)
    genotype_data[str(algo_prefix) + '_AD_TUMOR'] = str(nreads_ref_tumor) + ',' + str(nreads_alt_tumor)
    genotype_data[str(algo_prefix) + '_AD_CONTROL'] = str(nreads_ref_control) + ',' + str(nreads_alt_control)
    genotype_data[str(algo_prefix) + '_DP_TUMOR'] = str(nreads_tumor)
    genotype_data[str(algo_prefix) + '_DP_CONTROL'] = str(nreads_control) 

def merge_multiple_vcfs(tumor_sample_id, control_sample_id, merge_all_vcf, merge_somatic_vcf, logger, mutect_vcf, strelka_snv_vcf, strelka_indel_vcf, compress):

   all_calls = {}
   all_vcf_meta = {}
   all_vcf_meta_filter = {}

   if not mutect_vcf is None:
        merge_vcf_calls(mutect_vcf, tumor_sample_id, control_sample_id, all_calls, 'mutect2','snv')
        add_meta_info(all_vcf_meta, algorithm = 'mutect2')
   if not strelka_snv_vcf is None:
        merge_vcf_calls(strelka_snv_vcf, tumor_sample_id, control_sample_id, all_calls, 'strelka2', 'snv')
        add_meta_info(all_vcf_meta, algorithm = 'strelka2')
   if not strelka_indel_vcf is None:
        merge_vcf_calls(strelka_indel_vcf, tumor_sample_id, control_sample_id, all_calls, 'strelka2', 'indel')
        add_meta_info(all_vcf_meta, algorithm = 'strelka2')


   strelka_snvs_to_remove = []
   for chrom in chromosomes:
       if chrom in all_calls.keys():
           for varkey in all_calls[chrom].keys():
               varkey_elements = varkey.split('_')
               pos = varkey_elements[1]
               ref_allele = varkey_elements[2]
               alt_allele = varkey_elements[3]

               ## Check MuTect MNVs (multinucleotide variants/block substitutions) for support among consecutive Strelka (PASS) SNVs
               ##  - if all SNVs contributing to the MNVs are found as 'PASS' among Strelka calls, 
               ##    set 'Strelka' support for the candidate MNV and remove the individual Strelka SNVs
               if len(ref_allele) > 1 and len(alt_allele) > 1 and len(ref_allele) == len(alt_allele):
                   if 'PASS' in all_calls[chrom][varkey]['filter'].keys():
                        
                        len_allele = len(ref_allele)
                        ref_bases = list(ref_allele)
                        alt_bases = list(alt_allele)
                        i = 0

                        strelka_all_candidates = []
                        strelka_pass_candidates = []

                        while i < len_allele:
                            varkey_strelka = chrom + '_' + str(int(pos) + i) + '_' + str(ref_bases[i]) + '_' + str(alt_bases[i])

                            if varkey_strelka in all_calls[chrom].keys():
                                if 'strelka2' in all_calls[chrom][varkey_strelka]['callers_somatic'].keys():
                                    vfilters = ';'.join(all_calls[chrom][varkey_strelka]['filter'].keys())
                                    strelka_pass_candidates.append(varkey_strelka + ":" + str(vfilters))
                                else:
                                    vfilters = ';'.join(all_calls[chrom][varkey_strelka]['filter'].keys())
                                    strelka_all_candidates.append(varkey_strelka + ':' + str(vfilters))
                    
                            i = i + 1
                        
                        if len(strelka_pass_candidates) == len_allele:
                            all_calls[chrom][varkey]['callers_somatic']['strelka2'] = 1
                            all_calls[chrom][varkey]['info']['MNV_SUPPORT_STRELKA'] = "MNV_SUPPORT_STRELKA=" + ",".join(strelka_pass_candidates)
                            for c in strelka_pass_candidates:
                                strelka_snvs_to_remove.append(re.sub(r':PASS$', '', c))



   for chrom in chromosomes:
       if chrom in all_calls.keys():
            for var in strelka_snvs_to_remove:
                if var in all_calls[chrom].keys():
                    all_calls[chrom].pop(var)

   f = open(merge_all_vcf,'w')
   
   f.write("##fileformat=VCFv4.2\n")
   for k in sorted(all_vcf_meta.keys()):
       f.write(str(all_vcf_meta[k]) + '\n')
   for k in sorted(all_vcf_meta_filter.keys()):
       f.write(str(k) + '\n')
   f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
   
   f2 = open(merge_somatic_vcf,'w')
   f2.write("##fileformat=VCFv4.2\n")
   for k in sorted(all_vcf_meta.keys()):
       f2.write(str(all_vcf_meta[k]) + '\n')
   for k in sorted(all_vcf_meta_filter.keys()):
       f2.write(str(k) + '\n')
   f2.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

   n_somatic = 0
   n_somatic_strong = 0
   n_somatic_mutect = 0
   n_somatic_strelka = 0
   for chrom in chromosomes:
       if chrom in all_calls.keys():
            sorted_calls = {}
            sorted_calls_somatic = {}
            for varkey in all_calls[chrom]:
                fixed_str = str(chrom) + '\t' + str(all_calls[chrom][varkey]['pos']) + '\t.\t' + str(all_calls[chrom][varkey]['ref']) + '\t' + str(all_calls[chrom][varkey]['alt']) + '\t.\t'
                ref_allele = str(all_calls[chrom][varkey]['ref'])
                alt_allele = str(all_calls[chrom][varkey]['alt'])
                somatic = 0

                if "SOMATIC" in all_calls[chrom][varkey]['info'].keys():             
                    somatic = 1
                    n_somatic += 1
                    if len(all_calls[chrom][varkey]['callers_somatic'].keys()) == 2:
                        n_somatic_strong += 1                    
                    else:
                        if 'strelka2' in all_calls[chrom][varkey]['callers_somatic'].keys():
                            n_somatic_strelka += 1
                        else:
                            n_somatic_mutect += 1
   

                all_calls[chrom][varkey]['info']['VARIANT_CALLERS'] = 'VARIANT_CALLERS=' +  str(','.join(all_calls[chrom][varkey]['callers'].keys()))
                if somatic == 1:
                    all_calls[chrom][varkey]['info']['VARIANT_CALLERS'] = 'VARIANT_CALLERS=' +  str(','.join(all_calls[chrom][varkey]['callers_somatic'].keys()))
                if 'mutect2' in all_calls[chrom][varkey]['callers'].keys():
                    all_calls[chrom][varkey]['info']['TDP'] = 'TDP=' + str(all_calls[chrom][varkey]['info']['MTCT2_DP_TUMOR'].split('=')[1])
                    all_calls[chrom][varkey]['info']['TVAF'] = 'TVAF=' + str(all_calls[chrom][varkey]['info']['MTCT2_VAF_TUMOR'].split('=')[1])
                    all_calls[chrom][varkey]['info']['CDP'] = 'CDP=' + str(all_calls[chrom][varkey]['info']['MTCT2_DP_CONTROL'].split('=')[1])
                    all_calls[chrom][varkey]['info']['CVAF'] = 'CVAF=' + str(all_calls[chrom][varkey]['info']['MTCT2_VAF_CONTROL'].split('=')[1])
                else:
                    all_calls[chrom][varkey]['info']['TDP'] = 'TDP=' + str(all_calls[chrom][varkey]['info']['STKA2_DP_TUMOR'].split('=')[1])
                    all_calls[chrom][varkey]['info']['TVAF'] = 'TVAF=' + str(all_calls[chrom][varkey]['info']['STKA2_VAF_TUMOR'].split('=')[1])
                    all_calls[chrom][varkey]['info']['CDP'] = 'CDP=' + str(all_calls[chrom][varkey]['info']['STKA2_DP_CONTROL'].split('=')[1])
                    all_calls[chrom][varkey]['info']['CVAF'] = 'CVAF=' + str(all_calls[chrom][varkey]['info']['STKA2_VAF_CONTROL'].split('=')[1])

                all_calls[chrom][varkey]['info']['VT'] = 'VT=.'
                if(len(ref_allele) == 1 and len(ref_allele) == len(alt_allele)):
                    all_calls[chrom][varkey]['info']['VT'] = 'VT=snv'
                if(len(ref_allele) > 1 and len(alt_allele) > 1):
                    all_calls[chrom][varkey]['info']['VT'] = 'VT=blocksub'
                if(len(ref_allele) == 1 and len(alt_allele) > 1):
                    all_calls[chrom][varkey]['info']['VT'] = 'VT=INdel'
                if(len(ref_allele) > 1 and len(alt_allele) == 1):
                    all_calls[chrom][varkey]['info']['VT'] = 'VT=inDEL'
                info_str = ';'.join(map(str, sorted(all_calls[chrom][varkey]['info'].values())))
                filter_val = ';'.join(all_calls[chrom][varkey]['filter'].keys())
                record_line = str(fixed_str) + str(filter_val) + '\t' + str(info_str) + '\n'
                if not int(all_calls[chrom][varkey]['pos']) in sorted_calls.keys():
                    sorted_calls[int(all_calls[chrom][varkey]['pos'])] = []
                if somatic == 1 and not int(all_calls[chrom][varkey]['pos']) in sorted_calls_somatic.keys():
                    sorted_calls_somatic[int(all_calls[chrom][varkey]['pos'])] = []
                sorted_calls[int(all_calls[chrom][varkey]['pos'])].append(record_line)
                if somatic == 1:
                    record_line = str(fixed_str) + 'PASS' + '\t' + str(info_str) + '\n'
                    sorted_calls_somatic[int(all_calls[chrom][varkey]['pos'])].append(record_line)
            for pos in sorted(sorted_calls.keys()):
                for var_record in sorted_calls[pos]:
                    f.write(var_record)
            for pos in sorted(sorted_calls_somatic.keys()):
                for var_record in sorted_calls_somatic[pos]:
                    f2.write(var_record)

   logger.info('Total somatic calls: ' +  str(n_somatic) + ' / MuTect+Strelka: ' + str(n_somatic_strong) + ' / MuTect: ' + str(n_somatic_mutect) + ' / Strelka: ' + str(n_somatic_strelka))
   f.close()
   f2.close()

   if compress is True:
       os.system('bgzip -f ' + str(merge_all_vcf))
       os.system('tabix -p vcf ' + str(merge_all_vcf) + '.gz')
       os.system('bgzip -f ' + str(merge_somatic_vcf))
       os.system('tabix -p vcf ' + str(merge_somatic_vcf) + '.gz')

def merge_vcf_calls(vcffile, tumor_sample_id, control_sample_id, calls, algorithm = 'mutect2', vartype = 'snv'):

   if not vcffile is None:
       if os.path.exists(vcffile):
            for chrom in chromosomes:
                if not chrom in calls:
                    calls[chrom] = {}
        
            tumor_sample_index_column = -1
            control_sample_index_column = -1
            vcf_genotype_data_control = None
            vcf_genotype_data_tumor = None
            with gzip.open(vcffile,'rt') as fi:
                for line in fi:
                    if line.startswith('#CHROM'):
                        sample_column_indices = get_sample_column_index(line, tumor_sample_id, control_sample_id)
                        tumor_sample_index_column = sample_column_indices['tumor_sample_index_column']
                        control_sample_index_column = sample_column_indices['control_sample_index_column']
                    if not line.startswith('#'):
                        vcf_variant_record = line.rstrip().split('\t')
                        chromosome = str(re.sub(r'chr','',vcf_variant_record[0]))
                        pos = str(vcf_variant_record[1])
                        ref = str(vcf_variant_record[3])
                        alts = str(vcf_variant_record[4])
                        variant_filter = str(vcf_variant_record[6])

                        if ',' in alts and variant_filter == 'PASS':
                            print("WARNING: multi-allelic site with a PASS filter that is ignored")
                        ## Ignoring multiallelic sites for now
                        if alts == '.' or ',' in alts:
                            continue
                        #if alts == '.':
                            #continue
                        for alt in alts.split(','):
                            vcf_variant_key = chromosome + '_' + pos + '_' + ref + '_' + alt
                            genotype_data = init_genotype(algorithm = algorithm)
                            add_filter(algorithm, genotype_data, variant_filter)
                            vcf_info_data = get_info_data(vcf_variant_record[7].split(';'))
                            if len(vcf_variant_record) > 8 and tumor_sample_index_column != -1 and control_sample_index_column != -1:
                                format_tags = vcf_variant_record[8].split(':')
                                vcf_genotype_data_tumor = get_gt_data(vcf_variant_record[tumor_sample_index_column].split(':'), format_tags)
                                vcf_genotype_data_control = get_gt_data(vcf_variant_record[control_sample_index_column].split(':'), format_tags)
                    
                            populate_genotype_info_tags(genotype_data, ref, alt, vcf_info_data, vcf_genotype_data_tumor, vcf_genotype_data_control, algorithm = algorithm, vartype = vartype)

                            if not vcf_variant_key in calls[chromosome].keys():
                                calls[chromosome][vcf_variant_key] = {}
                                calls[chromosome][vcf_variant_key]['ref'] = ref
                                calls[chromosome][vcf_variant_key]['alt'] = alt
                                calls[chromosome][vcf_variant_key]['pos'] = pos
                                calls[chromosome][vcf_variant_key]['filter'] = {}
                                calls[chromosome][vcf_variant_key]['info'] = {}
                                calls[chromosome][vcf_variant_key]['callers'] = {}
                                calls[chromosome][vcf_variant_key]['callers_somatic'] = {}


                            for k in genotype_data.keys():
                                if not k == 'STKA2_FAIL_REASONS':
                                    calls[chromosome][vcf_variant_key]['info'][str(k)] = str(k) + '=' + str(genotype_data[k])
                                if 'FILTER'in k:
                                    calls[chromosome][vcf_variant_key]['info'][str(k)] = re.sub(r';',',',calls[chromosome][vcf_variant_key]['info'][str(k)])
                                    for f in genotype_data[k].split(';'):
                                        if not f in calls[chromosome][vcf_variant_key]['filter'].keys():
                                            calls[chromosome][vcf_variant_key]['filter'][f] = 1
                                        else:
                                            calls[chromosome][vcf_variant_key]['filter'][f] += 1
                                        calls[chromosome][vcf_variant_key]['callers'][algorithm] = 1                                       
                                        if f == 'PASS':
                                            calls[chromosome][vcf_variant_key]['callers_somatic'][algorithm] = 1
                                        
                                        if f == 'PASS':
                                            calls[chromosome][vcf_variant_key]['info']['SOMATIC'] = 'SOMATIC'
                                        else:
                                            if algorithm == 'strelka2':
                                                calls[chromosome][vcf_variant_key]['info']['STKA2_FAIL_REASONS'] = 'STKA2_FAIL_REASONS=' + str(f)

                                       
            fi.close()     
        
  
   return 0 

if __name__=="__main__": __main__()

         
