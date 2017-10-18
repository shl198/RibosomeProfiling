import re
import pandas as pd


class ncbi_gff(object):
    def __init__(self,df):
        self.df = df
        self.df.columns=['chr','source','feature','start','end','score','strand','frame','anno']
        self.df['start'] = self.df['start'] - 1
        self.df = self.df[self.df['feature'].values!='region']
        self.df = self.df.reset_index(drop=True)
        self.df['geneid'] = self.df['anno'].apply(lambda x: ncbi_gff.get_id(x,'GeneID:'))
        self.df['trid'] = self.df['anno'].apply(lambda x: ncbi_gff.get_id(x,'transcript_id='))
        self.df['prid'] = self.df['anno'].apply(lambda x: ncbi_gff.get_id(x,'protein_id='))
    @staticmethod
    def get_id(anno,feature):
        '''get id based on the feature provided'''
        try:
            gene_id = re.search('(?<={id}).+?(?=[;,]|$)'.format(id=feature),anno).group(0)
        except:
            gene_id = '-'
        return gene_id
    
    @staticmethod
    def get_tr_longest_intron(tr_df):
        '''get the longest intron in the transcript'''
        start = tr_df['start'].tolist()
        end = tr_df['end'].tolist()
        strand = tr_df['strand'].tolist()
        if len(start) == 1:
            return 0
        if strand[0] == '+':
            intron = max([abs(int(s)-int(e)) for s,e in zip(start[1:],end[:-1])])
        else:
            intron = max([abs(int(s)-int(e)) for s,e in zip(start[:-1],end[1:])])
        return intron
    
    def get_longest_intron(self):
        '''this is the longest intron across the whole genome'''
        df = self.df
        df = df[(~df['prid'].isnull()) & (df['feature'].values=='CDS')]
        df = df.reset_index(drop=True)
        df = df.groupby(['chr','prid']).apply(ncbi_gff.get_tr_longest_intron)
        return df#.max()
    
    def get_all_id(self,other_info=False):
        '''this function gets all ids in the gff file
            if other_info == True, output all columns of df
        '''
        df = self.df
        id_df = df[df['feature'].isin(['exon','CDS'])]
        id_df = id_df.reset_index(drop=True)
        
        id_df['sym'] = id_df['anno'].map(lambda x: ncbi_gff.get_id(x,'gene='))
        id_df['rna'] = id_df['anno'].map(lambda x: ncbi_gff.get_id(x,'Parent='))
        
        exn_df = id_df[id_df['feature'].values=='exon'][['geneid','sym','chr','rna','trid']].drop_duplicates()
        exn_df = exn_df.reset_index(drop=True)
        
        cds_df = id_df[id_df['feature'].values=='CDS'][['geneid','sym','chr','rna','prid']].drop_duplicates()
        cds_df = cds_df.reset_index(drop=True)
        
        merge_df = pd.merge(exn_df,cds_df,how='outer',on=['geneid','sym','chr','rna'])
        merge_df.columns = ['GeneID','GeneSymbol','Chrom','TrID','TrAccess','PrAccess']
        merge_df.fillna('-',inplace=True)
        merge_df.fillna('-',inplace=True)
#         merge_df = merge_df[(merge_df['TrAccess'].values != '-') | (merge_df['PrAccess'].values != '-')]
        merge_df = merge_df.sort_values(['GeneID'])
        if other_info == True:
            return merge_df
        else:
            return merge_df[['GeneID','GeneSymbol','Chrom','TrAccess','PrAccess','TrID']]
    
    def get_gene_seq(self,ref_dic,gid,id_type='tr'):
        '''this function gets seqeunce of a transcript or protein
        * ref_dic: handle reads using SeqIO.index(fa,'fasta')
        * gid: gene id
        '''
        df = self.df
        if id_type == 'tr':
            feature = 'exon'
            id_t = 'trid'
        elif id_type == 'pr':
            feature = 'CDS'
            id_t = 'prid'
        region_df = df[(df['feature'].values==feature) & (df[id_t].values==gid)]
        # get sequence
        scaff = region_df['chr'].tolist()[0]
        scaff_seq = ref_dic[scaff].seq
        strand = region_df['strand'].tolist()[0]
        
        g_seq = ''
        for s,e in zip(region_df['start'],region_df['end']):
            g_seq += scaff_seq[int(s):int(e)]
        # consider strand
        if strand == '-':
            g_seq = g_seq.reverse_complement()
        
        if id_type == 'pr':
            g_seq = g_seq.translate()
        return g_seq
    
    def get_feature_pos_df(self,feature):
        '''this gets either cds or exon position dataframes'''
        df = self.df
        df = df[df['feature'].values==feature]
        df = df.reset_index(drop=True)

        if feature=='CDS':
            res_df = df[['chr','start','end','geneid','prid','strand']]
            res_df = res_df.dropna()
            res_df.columns = ['Chr', 'cds_s','cds_e','GeneID','PrAccess','Strand']

        elif feature == 'exon':
            res_df = df[['chr','start','end','geneid','trid','strand']]
            res_df = res_df.dropna()
            res_df.columns = ['Chr', 'exn_s','exn_e','GeneID','TrAccess','Strand']
        return res_df

    def multi_chr_protein(self):
        proteins = []
        df = self.df
        chr_gene_df = df[['chr','prid']].drop_duplicates()
        gene_chr_dic = {k:list(v) for k,v in chr_gene_df.groupby('prid')['chr']}
        for key in gene_chr_dic:
            if len(set(gene_chr_dic[key])) > 1:
                proteins.append(key)
        return proteins



# gff_fn = '/data/shangzhong/RibosomeProfiling/ritux/reference/combine.gff'
# df = pd.read_csv(gff_fn,sep='\t',header=None,comment='#')
# obj = ncbi_gff(df)
# all_id_df = obj.get_all_id()
# all_id_df.to_csv('/data/shangzhong/RibosomeProfiling/ritux/reference/combine_id.txt',sep='\t',index=False)
# 
# ref_dic = SeqIO.index('/data/genome/hamster/picr/picr.fa','fasta')
# res = obj.get_gene_seq(ref_dic,'NM_001246795',id_type='tr')
# print res