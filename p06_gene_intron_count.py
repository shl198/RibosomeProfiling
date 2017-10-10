import os,argparse,glob
import pandas as pd
from natsort import natsorted
from multiprocessing import Pool



parser = argparse.ArgumentParser(description='Get read count at cds, utr5 and utr3 for each protein')
parser.add_argument('-g','--gff_path',action='store',dest='gff_path',help='path to all files generated by p01_prepare_annotation.py')
parser.add_argument('-c','--cov_path',action='store',dest='cov_path',help='coverage path that stores all the coverage files')
parser.add_argument('-o','--offset',action='store',dest='offset',help='full path to offset file')
parser.add_argument('-t','--thread',action='store',dest='thread',type=int,help='thread to use',default=1)
parser.add_argument('-e','--end',action='store',dest='end',help='which end of read to use',default=5)
args = parser.parse_args()


db_path = args.gff_path
if db_path.endswith('/'): db_path = db_path[:-1]
cov_path = args.cov_path
if cov_path.endswith('/'): cov_path = cov_path[:-1]
offset_fn = args.offset
end = args.end
thread = args.thread





def get_pos_dic(bed):
    '''this function builds dictioanry for protein or rna {id:['+/-',(s1,e1),(s2,e2)]}'''
    dic = {}
    with open(bed) as f:
        for line in f:
            if line.startswith('Chr'): continue
            item = line.strip().split('\t')
            if item[4] in dic:
                dic[item[4]].append((item[1],item[2]))
            else:
                dic[item[4]] = [item[-1],(item[1],item[2])]
    return dic

def get_pos(dic,access):
    '''dic format is {access:['+/-',(s1,e1),(s2,e2)]}'''
    if access not in dic:
        pos = []
    else:
        strd = dic[access][0]
        pos = []
        for p in dic[access][1:]:
            pos.extend(range(int(p[0])+1,int(p[1])+1))
        if strd == '-':
            pos = pos[::-1]
    return pos


def cov5_3_dic(covFile,m_lens):
    '''prepare two dictionaries for mapping of 5end and 3end of the reads.
    format {chr:{pos+/-:count}}'''
    cov_5dic = {}
    cov_3dic = {}
    with open(covFile) as cov:
        for line in cov:
            item = line.strip().split('\t')
            if int(item[-1]) not in m_lens: continue
            count = int(item[0])
            chrom = item[1]
            end5  = item[2]
            end3  = item[3]
            strd  = item[4]
            if chrom in cov_5dic:
                if end5+strd in cov_5dic[chrom]:
                    cov_5dic[chrom][end5+strd] += count
                else:
                    cov_5dic[chrom][end5+strd] = count
                if end3+strd in cov_3dic[chrom]:
                    cov_3dic[chrom][end3+strd] += count
                else:
                    cov_3dic[chrom][end3+strd] = count
            else:
                cov_5dic[chrom] = {}
                cov_3dic[chrom] = {}
    return cov_5dic,cov_3dic


def get_pos_cov(dic5,dic3,chrom,pos,end,strand='no'):
    '''if does not proved strand information it use the sequence position to decide'''
    end = str(end)
    if strand == '+' or strand == '-':
        strd = strand
    else:
        if pos[0]<pos[1]:
            strd = '+'
        else:
            strd = '-'
    pos_cov = []
    for p in pos:
        try:
            if end == '5' and strd == '+':
                pos_cov.append(dic5[chrom][str(p)+strd])
            elif end == '5' and strd == '-':
                pos_cov.append(dic3[chrom][str(p)+strd])
            elif end == '3' and strd == '+':
                pos_cov.append(dic3[chrom][str(p)+strd])
            elif end == '3' and strd == '-':
                pos_cov.append(dic5[chrom][str(p)+strd])
        except:
            pos_cov.append(0)
    return pos_cov

def single_gene_count(gene,ge_pr_dic,utr_df,dic5,dic3,exn_pos_dic,cds_pos_dic,offset,end,intron=True):
    '''
    * dic5,dic3: coverage dictionary using 5/3 end of reads
    * end: which end of reads to consider to count reads
    '''
    gene_pos_dic = {} # store all position for a gene
    if str(end) == '3':
        offset = - offset
    gene_pos = []
    prs = ge_pr_dic[gene]
    for pr in prs:
        tr = utr_df.loc[pr,'TrAccess']
        utr5  = int(utr_df.loc[pr,'utr5_len'])
        utr3  = int(utr_df.loc[pr,'utr3_len'])
        chrom = utr_df.loc[pr,'Chrom']
        strand = utr_df.loc[pr,'strand']
        # get the rna position
        if tr in exn_pos_dic:
            tr_pos = get_pos(exn_pos_dic,tr)
        else:
            tr_pos = get_pos(cds_pos_dic,pr)
        s = tr_pos[0];e = tr_pos[-1]
        shift = abs(offset)
        if tr_pos[0] < tr_pos[1]:
            new_pos = list(range(s-shift,s)) + tr_pos + list(range(e+1,e+shift+1))
        else:
            new_pos = list(range(s+shift,s,-1)) + tr_pos + list(range(e-1,e-1-shift,-1))
        cds_pos = new_pos[shift+utr5-offset:len(new_pos)-shift-utr3-offset]
        # add position to chromosome
        if chrom not in gene_pos_dic:
            gene_pos_dic[chrom] = [[],[]] # one for positive, one for negative
        if strand == '+':
            gene_pos_dic[chrom][0].extend(cds_pos)
        else:
            gene_pos_dic[chrom][1].extend(cds_pos)
    # get coverage
    gene_count = 0
    intron_count = 0
    for chrom in gene_pos_dic:
        for pos in gene_pos_dic[chrom]:
            if pos != []:
                pos = set(pos)
                if intron == True:
                    intron_pos = [p for p in list(range(min(pos),max(pos)+1)) if p not in pos]
                    if intron_pos != []:
                    	intron_count += sum(get_pos_cov(dic5,dic3,chrom,intron_pos,end,strand))
                gene_count += sum(get_pos_cov(dic5,dic3,chrom,list(pos),end,strand))
    return [gene_count,intron_count]
    



def get_gene_intron_count(gene_count_path,covFile,cdsFile,exnFile,utr_fn,offset_fn,end):
    '''get total count of each gene
    * gene_count_path: output path that stores all the output files'''
    # prepare position dictionary
    cds_pos_dic = get_pos_dic(cdsFile)
    exn_pos_dic = get_pos_dic(exnFile)
    # read utr data, get proteins and build dictionaries
    utr_df = pd.read_csv(utr_fn,sep='\t',header=0)
    utr_df.index = utr_df['PrAccess']
    utr_df['utr5_len'] = utr_df['utr5_len'].astype('int')
    utr_df['utr3_len'] = utr_df['utr3_len'].astype('int')
    ge_pr_dic = {k:list(v) for k,v in utr_df.groupby('GeneID')['PrAccess']}
    # read offset file
    offset_df = pd.read_csv(offset_fn,sep='\t',header=0)
    offset_dic = {k:list(v) for k,v in offset_df.groupby('offset')['len']}
    # count reads
    intron = True
    gene_count = {}
    intron_count = {}

    for offset in offset_dic:
        # prepare coverage dictionary
        dic5,dic3 = cov5_3_dic(covFile,offset_dic[offset])
        for gene in ge_pr_dic:
            if gene not in gene_count:
                gene_count[gene]   = {}
                intron_count[gene] = {}
            count = single_gene_count(gene,ge_pr_dic,utr_df,dic5,dic3,exn_pos_dic,cds_pos_dic,offset,str(end),True)
            gene_count[gene][offset] = count[0]
            if intron ==True:
                intron_count[gene][offset] = count[1]
    outFile = gene_count_path + '/' + os.path.basename(covFile)
    # output to file
    with open(outFile,'w') as f:
        f.write('\t'.join(['geneid','gene_count','intron_count'])+'\n')
        for gene in gene_count:
            g_c = sum(gene_count[gene].itervalues())  # in python3 use .values instead of .itervalues
            i_c = sum(intron_count[gene].itervalues())
            f.write('\t'.join([gene,str(g_c),str(i_c)]) + '\n')        


import time
start = time.time()


exnFile = db_path + '/01_pr_rna.bed'
cdsFile = db_path + '/01_pr_cds.bed'
utr_fn = db_path + '/03_utr_len.txt'
# offset_fn = db_path + '/p_offset.txt'

path = os.path.dirname(db_path)
cov_path = path + '/02_cov'
covFiles = natsorted(glob.glob(cov_path+'/*.txt'))
gene_count_path = path + '/06_gene_intron_count'
if not os.path.exists(gene_count_path): os.mkdir(gene_count_path)


# p = Pool(processes=thread)
# for covFile in covFiles:
# 	p.apply_async(get_gene_intron_count,args=(gene_count_path,covFile,cdsFile,exnFile,utr_fn,offset_fn,end,))
# p.close()
# p.join()

get_gene_intron_count(gene_count_path,covFiles[0],cdsFile,exnFile,utr_fn,offset_fn,end)
print('get gene intron count succeed')
print('total run time is: '+str((time.time()-start)/60) +' minutes')
