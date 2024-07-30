import os
import argparse
import numpy as np
import multiprocessing


parser=argparse.ArgumentParser()
parser.add_argument("-kmer_size",help="Size of kmer. Increasing this number quickly increases memory consumption.",required=True,type=int,metavar='[int]')
parser.add_argument("-fasta_folder",help="Path to folder containing (only) the fasta files to reduce to kmer frequencies.",required=True,metavar='[str]')
parser.add_argument("-output_folder",help="Path to output folder to store kmer frequencies and index file. Created if missing",required=True,metavar='[str]')
parser.add_argument("-p",help="Number of processes for parallelization. Setting a number higher than the number of fasta files is useless. Default: 1",type=int,default=1,metavar='[int]')
parser.add_argument("-output_type",help="Type of output. single_files: one numpy array for each fasta file. matrix: one numpy matrix where each row stores the kmer frequencies of the corresponding fasta file as indicated in the output index. contig_matrix: one numpy matrix where each row stores the kmer frequencies of the corresponding contig as indicated in the output index. All contigs from all fasta files are stored in this matrix. Default: single_files",choices=['single_files','fasta_matrix','contig_matrix'],default='single_files')
parser.add_argument("--get_kmers",help="Setting this flag will print to stdout a comma-separated list of kmers in the same order of the frequencies in the output numpy arrays. The program will terminate before performing other operations.",action='store_true')
args=parser.parse_args()

def chunks(l,nproc):
    right_div=len(l)//nproc
    nmore=len(l)%nproc
    return [l[i*(right_div+1):i*(right_div+1)+right_div+1] for i in range(nmore)]+[l[nmore*(right_div+1)+i*right_div:nmore*(right_div+1)+i*right_div+right_div] for i in range(nproc-nmore) if nmore<len(l)]

def invcomp(seq):
    toret=''
    for j in range(len(seq)-1,-1,-1):
        i=seq[j]
        if i=='A':
            toret+='T'
        elif i=='C':
            toret+='G'
        elif i=='G':
            toret+='C'
        elif i=='T':
            toret+='A'
    return toret

def all_kmers(ksize,kmer,diz): 
    if len(kmer)<ksize:
        for i in 'ACGT':
            all_kmers(ksize,kmer+i,diz)
    else:
        diz[kmer]=0

def sub_allkmers(d):
    newd={}
    for i in d:
        ic=invcomp(i)
        if i not in newd and ic not in newd:
            newd[i]=0
    return newd
                
def contigsToFreqs(ksize,contigs):
    kcount={}
    allowed=['A','C','G','T']
    all_kmers(ksize,'',kcount)
    kcount=sub_allkmers(kcount)
    for contig in contigs:
        for i in range(len(contig)-ksize+1):
            kmer=contig[i:i+ksize]
            if not all([base in allowed for base in kmer]): # skip kmers with N or masked sequences
                continue
            if kmer in kcount: kcount[kmer]+=1
            else: kcount[invcomp(kmer)]+=1
    totkmers=sum(kcount.values())
    return np.array([kcount[x]/totkmers for x in kcount])

def contigToFreqs(ksize,contig):
    kcount={}
    allowed=['A','C','G','T']
    all_kmers(ksize,'',kcount)
    kcount=sub_allkmers(kcount)
    for i in range(len(contig)-ksize+1):
        kmer=contig[i:i+ksize]
        if not all([base in allowed for base in kmer]): # skip kmers with N or masked sequences
            continue
        if kmer in kcount: kcount[kmer]+=1
        else: kcount[invcomp(kmer)]+=1
    totkmers=sum(kcount.values())
    return np.array([kcount[x]/totkmers for x in kcount])

def main(args):
    fasta_files,ksize,outfold,infold=args
    toret=[]
    for fasta in fasta_files:
        contigs=[]
        f=open(infold+'/'+fasta)
        contig=[]
        for line in f:
            if line[0]=='>' and contig:
                contigs.append(''.join(contig))
                contig=[]
            else:
                contig.append(line.strip())
        contigs.append(''.join(contig))
        f.close()
        freqs=contigsToFreqs(ksize,contigs)
        toret.append(freqs)
    return toret

def main_contig(args):
    fasta_files,ksize,outfold,infold=args
    toret=[]
    names=[]
    for fasta in fasta_files:
        f=open(infold+'/'+fasta)
        contig=[]
        cname=f.readline()[1:].strip().split(' ')[0]
        for line in f:
            if line[0]=='>':
                freqs=contigToFreqs(ksize,''.join(contig))
                toret.append(freqs)
                names.append(cname)
                cname=line[1:].strip().split(' ')[0]
                contig=[]
            else:
                contig.append(line.strip())
        freqs=contigToFreqs(ksize,''.join(contig))
        toret.append(freqs)
        names.append(cname)
        f.close()
    return toret,names
        

if __name__=='__main__':

    if args.get_kmers:
        kcount={}
        all_kmers(args.kmer_size,'',kcount)
        kcount=sub_allkmers(kcount)
        print(','.join(x for x in kcount))
        exit()

    # check output folder
    if not os.path.isdir(args.output_folder):
        if os.path.exists(args.output_folder): print('Error:',args.output_folder,'exists and is not a directory.\nAborting...'); exit(1)
        else: os.mkdir(args.output_folder)
    # prepare fasta files for parallelization
    fasta_files=os.listdir(args.fasta_folder)
    chunked_fastas=chunks(fasta_files,args.p)
    # parallel reduction to kmer frequencies
    print('Starting to reduce to kmer frequencies. This will take a while...')
    if args.output_type!='contig_matrix':
        with multiprocessing.Pool(args.p) as pool:
            all_freqs=pool.map(main,[(fastas,args.kmer_size,args.output_folder,args.fasta_folder) for fastas in chunked_fastas])
        print('Saving output array(s)...')
        toconcat,index=[],[]
        for x,y in zip([x for y in all_freqs for x in y],[x for y in chunked_fastas for x in y]):
            if args.output_type=='single_files': np.save(args.output_folder+'/'+y,x)
            elif args.output_type=='fasta_matrix': toconcat.append(np.expand_dims(x,0)); index.append(y)
        if args.output_type=='fasta_matrix': 
            np.save(args.output_folder+'/fasta_matrix',np.concatenate(toconcat,0))
            new=open(args.output_folder+'/fasta_matrix_index.txt','w'); print('\n'.join(index),file=new); new.close()
    else:
        with multiprocessing.Pool(args.p) as pool:
            res=pool.map(main_contig,[(fastas,args.kmer_size,args.output_folder,args.fasta_folder) for fastas in chunked_fastas])
        print('Saving output array(s)...')
        toconcat,index=[],[]
        for freqlist,namelist in res:
            for x,y in zip(freqlist,namelist):
                toconcat.append(np.expand_dims(x,0)); index.append(y)
        np.save(args.output_folder+'/contig_matrix',np.concatenate(toconcat,0))
        new=open(args.output_folder+'/contig_matrix_index.txt','w'); print('\n'.join(index),file=new); new.close()

    print('\nAll done! Output files are in '+args.output_folder+'.')