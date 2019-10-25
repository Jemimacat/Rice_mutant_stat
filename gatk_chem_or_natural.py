'''
This Script is used for variance statistics of freebayes VCFs
author: Li WENHUI, jemimalwh@gmail.com, liwenhui@genomics.cn
date: 2019-10-25
'''
import os, sys
import re

head = ['MTYPE','GROUP','CHANGE_BYTE','VARIANCE_NUM','VARIANCE_INFO']

def make_dictory(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(output_dir + '/data'):
        os.makedirs(output_dir + '/data')
    if not os.path.exists(output_dir + '/figure'):
        os.makedirs(output_dir + '/figure')

def ls_stat_files(input_dir):
    snp_files = []
    indel_file = []
    for f in os.listdir(input_dir):
        if 'snp' in f:
            snp_files.append(input_dir + '/' + f)
        elif 'indel' in f:
            indel_file.append(input_dir + '/' + f)
    return snp_files, indel_file

def open_f_handles(files):
    f_handles = []
    for file in files:
        f = open(file,'r')
        f_handles.append(f)   
    return f_handles

def close_f_handles(f_handles):
    for f in f_handles:
        f.close()

def variance_statistics(files,output_dir,vtype='snp',chem_th=50,nat_th=500):
    stat = {}
    stat['chem'] = {}
    stat['nature'] = {}
    ## open files
    f_handles = open_f_handles(files)
    for handle in f_handles:
        handle.readline()
        while(1):
            line = handle.readline().replace('\n','')
            if not line:
                break
            items = line.split('\t')
            mutant_num = int(items[3]) + int(items[4])
            chrom, pos, ref, alt = items[2].split('_')
            change_byte = '_'.join([ref,alt])
            if mutant_num <= 50:
                if not change_byte in stat['chem'].keys():
                    stat['chem'][change_byte] = {}
                stat['chem'][change_byte][items[2]] = mutant_num
            elif mutant_num <= 500:
                if not change_byte in stat['nature'].keys():
                     stat['nature'][change_byte] = {}
                stat['nature'][change_byte][items[2]] = mutant_num
    ## close files
    close_f_handles(f_handles)    
    ## output stat results
    for group in ('chem','nature'):
        f_out_file = output_dir + '/data/' + vtype + '.' + group + '.stat.tsv'
        f_out_handle = open(f_out_file,'w')
        f_out_handle.writelines('\t'.join(head)+'\n')
        for change_byte in stat[group].keys():
            sites = sorted(stat[group][change_byte].keys())
            num = len(sites)
            sites_info = ';'.join(sites)
            f_out_handle.writelines('\t'.join([vtype,group,change_byte,str(num),sites_info]) + '\n')
        f_out_handle.close()
        ## clean memory
        stat[group] = {}



if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    ## make directory
    make_dictory(output_dir)
    ## get snp & indel files
    snp_files, indel_files = ls_stat_files(input_dir)
    ## snp data stat
    variance_statistics(snp_files,output_dir,vtype='snp',chem_th=50,nat_th=500)
    ## indel data stat
    variance_statistics(indel_files,output_dir,vtype='indel',chem_th=50,nat_th=500)

