'''
This Script is used for variance statistics of freebayes VCFs
author: Li WENHUI, jemimalwh@gmail.com, liwenhui@genomics.cn
date: 2019-10-25
'''
import os, sys
import gzip
import re

## Reference information of Rice
chromosomes = ['chr01','chr02','chr03','chr04','chr05','chr06','chr07','chr08','chr09','chr10','chr11','chr12']
coordinates = [43270923,35937250,36413819,35502694,29958434,31248787,29697621,28443022,23012720,23207287,29021106,27531856]
reference = dict(zip(chromosomes,coordinates))
bin_width = 1000000
stat_results = {}

def make_bin_width(reference,bin_width=1000000):
    bin_widths = {}
    for chromosome in reference.keys():
        coordinate = reference[chromosome]
        bin_widths[chromosome] = []
        times = int(coordinate/bin_width)
        for i in range(0,times+1):
            start = i*bin_width + 1
            end = start + bin_width
            if end >= coordinate+1:
                end = coordinate+1
            bin_widths[chromosome].append([start,end])
    return bin_widths

bin_widths = make_bin_width(reference,bin_width=100000)

def print_head(outdir):
    head_file = outdir + '/header'
    out_head = ['TYPE','CHR','SITE','HOMOREF','HOMOALT','HYBRID','Homoref_info','Homoalt_info','Hybrid_info']
    f_head = open(head_file,'w')
    f_head.writelines('\t'.join(out_head) + '\n')
    f_head.close()

def make_dictory(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)   

def traverse_dir(input_dir):
    files = []
    for f in os.listdir(input_dir):
        files.append(input_dir + '/' + f)
    return files

def gzip_open_files(files):
    f_handles = []
    for file in files:
        #f = gzip.open(file,'rb')
        f = open(file,'r')
        f_handles.append(f)   
    return f_handles

def close_f_handles(f_handles):
    for f in f_handles:
        f.close()

def read_files_in_bathes(f_handles,outdir):

    def header_or_record(line):
        # header(1) or record(2)
        flag = 0
        h_pattern = re.compile(r'#CHROM')
        r_pattern = re.compile(r'chr')
        m1 = re.match(h_pattern,line)
        if (m1):
            flag = 1
        m2 = re.match(r_pattern,line)
        if (m2):
            flag = 2
        return flag

    def get_headers():
        headers = []
        for handle in f_handles:
            while(1):
                #line = handle.readline().decode().replace('\n','')
                line = handle.readline().replace('\n','')
                flag = header_or_record(line)
                if flag == 1:
                    headers.append(line.split('\t'))
                    break
        return headers

    def site_pass_or_not(info):
        pass_flag = 1
        return pass_flag

    def get_info(ref,alt,info):
        ## for freebayes vcf
        alt_ad_dep = []
        mutant = 'hybrid'
        mtype = 'indel'
        gt,dep,ad = info.split(':')[0:3]
        ads = ad.split(',')
        ads.pop(0)
        alts = alt.split(',')

        if '1/1' in gt:
            mutant = 'homoalt'
        elif '0/0' in gt:
            mutant = 'homoref'
        
        for i in range(0,len(ads)):
            this_ad = ads[i]
            this_alt = alts[i]
            if len(ref)==1 and len(this_alt)==1:
                mtype = 'snp'
            alt_ad_dep.append([mtype,mutant,this_alt,dep,this_ad])
        return alt_ad_dep

    def travers_handles(chromosome,start,end,bin_width):
        tag = '-'.join([str(start),str(end-1)])
        if not chromosome in data.keys():
            data[chromosome] = {}
        if not tag in data[chromosome].keys():
            data[chromosome][tag] = {}
        for i in range(0,len(f_handles)):
            handle = f_handles[i]
            header = headers[i]
            while(1):
                #line = handle.readline().decode().replace('\n','')
                line = handle.readline().replace('\n','')
                if not line:
                    break
                this_flag = header_or_record(line)
                if this_flag == 2:
                    temp = line.split('\t')
                    if site_pass_or_not(temp[7]):
                        if (temp[0]==chromosome) and (int(temp[1])>=start) and (int(temp[1])<end):                   
                            for j in range(9,len(temp)):
                                if temp[j] == '.':
                                    continue
                                samp_id = header[j]
                                alt_ad_dep = get_info(temp[3],temp[4],temp[j])
                                for mtype,mutant,alt,dep,ad in alt_ad_dep:
                                    if not mtype in data[temp[0]][tag].keys():
                                        data[temp[0]][tag][mtype] = {}
                                    if not '_'.join([temp[0],temp[1],temp[3],alt]) in data[temp[0]][tag][mtype].keys():
                                        data[temp[0]][tag][mtype]['_'.join([temp[0],temp[1],temp[3],alt])] = {}
                                    if not mutant in data[temp[0]][tag][mtype]['_'.join([temp[0],temp[1],temp[3],alt])].keys():
                                        data[temp[0]][tag][mtype]['_'.join([temp[0],temp[1],temp[3],alt])][mutant] = {}
                                    data[temp[0]][tag][mtype]['_'.join([temp[0],temp[1],temp[3],alt])][mutant][':'.join([samp_id,','.join([str(dep),str(ad)])])] = 1
                        else:
                            if not temp[0] in data.keys():
                                data[temp[0]] = {}
                            new_start = int(int(temp[1])/bin_width)*bin_width + 1
                            new_end = new_start + bin_width
                            new_tag = '-'.join([str(new_start),str(new_end-1)])
                            if not new_tag in data[temp[0]].keys():
                                data[temp[0]][new_tag] = {}
                            for j in range(9,len(temp)):
                                if temp[j] == '.':
                                    continue
                                samp_id = header[j]
                                alt_ad_dep = get_info(temp[3],temp[4],temp[j])
                                for mtype,mutant,alt,dep,ad in alt_ad_dep:
                                    if not mtype in data[temp[0]][new_tag].keys():
                                        data[temp[0]][new_tag][mtype] = {}
                                    if not '_'.join([temp[0],temp[1],temp[3],alt]) in data[temp[0]][new_tag][mtype].keys():
                                        data[temp[0]][new_tag][mtype]['_'.join([temp[0],temp[1],temp[3],alt])] = {}
                                    if not mutant in data[temp[0]][new_tag][mtype]['_'.join([temp[0],temp[1],temp[3],alt])].keys():
                                        data[temp[0]][new_tag][mtype]['_'.join([temp[0],temp[1],temp[3],alt])][mutant] = {}
                                    data[temp[0]][new_tag][mtype]['_'.join([temp[0],temp[1],temp[3],alt])][mutant][':'.join([samp_id,','.join([str(dep),str(ad)])])] = 1
                            break
        # stat variance of this Chromosome and tag, delele it, and keep the new tag
        stat_variance_in_bath(chromosome,tag)
        del data[chromosome][tag]

    def stat_variance_in_bath(chromosome,tag):
        batch_data = data[chromosome][tag]
        for mtype in batch_data.keys():
            fname = outdir + '/' + chromosome + '.'+ tag +'.' + mtype + '.stat.tsv'
            f_out = open(fname,'w')
            for site in batch_data[mtype].keys():
                homoref_info, homoalt_info, hybrid_info = ('-','-','-')
                homoref_num, homoalt_num, hybrid_num = (0,0,0)
                if 'homoref' in batch_data[mtype][site].keys():
                    homoref = sorted(batch_data[mtype][site]['homoref'].keys())
                    homoref_info = ';'.join(homoref)
                    homoref_num = len(homoref)
                if 'homoalt' in batch_data[mtype][site].keys():
                    homoalt = sorted(batch_data[mtype][site]['homoalt'].keys())
                    homoalt_info = ';'.join(homoalt)
                    homoalt_num = len(homoalt)
                if 'hybrid' in batch_data[mtype][site].keys():
                    hybrid = sorted(batch_data[mtype][site]['hybrid'].keys())
                    hybrid_info = ';'.join(hybrid)
                    hybrid_num = len(hybrid)
                f_out.writelines('\t'.join([mtype,chromosome,site,str(homoref_num),str(homoalt_num),str(hybrid_num),homoref_info,homoalt_info,hybrid_info]) + '\n')
            f_out.close()
            save_file_names(chromosome,mtype,fname)
    
    def save_file_names(chromosome,mtype,fname):
        if not chromosome in stat_results.keys():
            stat_results[chromosome] = {}
        if not mtype in stat_results[chromosome].keys():
            stat_results[chromosome][mtype] = []
        stat_results[chromosome][mtype].append(fname)
  
    def final_stat_variance():
        for chrom in list(data.keys()):
            for tag in list(data[chrom].keys()):
                for mtype in list(data[chrom][tag].keys()):
                    fname = outdir + '/' + chrom + '.' + tag + '.' + mtype + '.stat.tsv'
                    f_out = open(fname,'w')
                    for site in list(data[chrom][tag][mtype].keys()):
                        homoref_info, homoalt_info, hybrid_info = ('-','-','-')
                        homoref_num, homoalt_num, hybrid_num = (0,0,0)
                        if 'homoref' in data[chrom][tag][mtype][site].keys():
                            homoref = sorted(data[chrom][tag][mtype][site]['homoref'].keys())
                            homoref_info = ';'.join(homoref)
                            homoref_num = len(homoref)
                        if 'homoalt' in data[chrom][tag][mtype][site].keys():
                            homoalt = sorted(data[chrom][tag][mtype][site]['homoalt'].keys())
                            homoalt_info = ';'.join(homoalt)
                            homoalt_num = len(homoalt)
                        if 'hybrid' in data[chrom][tag][mtype][site].keys():
                            hybrid = sorted(data[chrom][tag][mtype][site]['hybrid'].keys())
                            hybrid_info = ';'.join(hybrid)
                            hybrid_num = len(hybrid)
                        f_out.writelines('\t'.join([mtype,chrom,site,str(homoref_num),str(homoalt_num),str(hybrid_num),homoref_info,homoalt_info,hybrid_info]) + '\n')
                    f_out.close()
                    save_file_names(chrom,mtype,fname)
                del data[chrom][tag]
    
    # Skip the "##" lines and keep the headers
    headers = get_headers()
    # Dictionary to archive data
    data = {}
    for chromosome in chromosomes:
        coordinations = bin_widths[chromosome]
        for start,end in coordinations:
            travers_handles(chromosome,start,end,bin_width)
    # Final clean archive
    final_stat_variance()
 
def combine_results(outdir):
    head_file = outdir + '/header'
    for chromosome in stat_results.keys():
        for mtype in stat_results[chromosome].keys():
            merge_fname = outdir + '/' + chromosome + '.' + mtype + '.tsv'
            fnames = ' '.join(stat_results[chromosome][mtype])
            os.system('cat ' + head_file + ' ' + fnames + ' > ' + merge_fname)
            os.system('rm -f ' + fnames)
    os.system('rm -f ' + head_file)


if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    ## make directory
    make_dictory(output_dir)
    ## print head to file
    print_head(output_dir)
    ## list files
    files = traverse_dir(input_dir)
    ## open all the files
    f_handles = gzip_open_files(files)
    ## read files and do the statistics in bathes
    read_files_in_bathes(f_handles,output_dir)
    ## merger results
    combine_results(output_dir)
    ## close all file handles
    close_f_handles(f_handles)
    
