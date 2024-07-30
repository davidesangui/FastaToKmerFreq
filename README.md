# FastaToKmerFreq
Simple script to reduce DNA sequences (in Fasta format) to K-mer frequencies.
## Installation
To use FastaToKmerFreq, simply download the file FastaToKmerFreq.py from this repository and run it with Python3.
The script requires numpy, which can be obtained via pip:
```
pip install numpy
```
## How it works
FastaToKmerFreq takes as input a folder of fasta files. All the contigs in each fasta file are reduced to a single set of k-mer frequencies representing the entire sequence stored in the the fasta file. The k-mer size is user-defined.  
In theory, there are 4^kmer_size possible kmers. But when considering forward and reverse k-mers and palindromic k-mers, there is lower number of distinct k-mers (e.g with k=4, there are 136 distinct kmers; for more details: https://stackoverflow.com/questions/40952719/algorithm-to-collapse-forward-and-reverse-complement-of-a-dna-sequence-in-python).  
FastaToKmerFreq reduces a sequence to a set of distinct k-mer frequencies by counting the k-mers in both the forward and reverse sequence.  
All k-mers containing characters other than {'A','C','G','T'} are ignored.  
The output is a folder of Numpy arrays generated with the function numpy.save(). Each output file is a Numpy array storing the k-mer frequencies of the corresponding sequence in fasta format.  
When multiple cores are available, different Fasta files can be reduced in parallel.
## Usage
FastaToKmerFreq.py requires three mandatory arguments:
- ```-kmer_size``` Size of k-mers. Since all possible k-mers are stored in memory at once, increasing the k-mer size quickly increases memory consumption. A rough estimation of memory consumption can be performed with the following formula: 8 * numCores * (4^kmer_size)/2
- ```-fasta_folder``` Path to folder containing (only) the fasta files to reduce to kmer frequencies.
- ```-output_folder``` Path to output folder to store kmer frequencies and index file. Created if missing.

**Example:**  

```python3 FastaToKmerFreq.py -kmer_size 3 -fasta_folder .\genomes\ -output_folder example_output -p 3``` 

```-p``` defines the number of processes for parallelization. Setting a number greater than the number of genomes is useless.  

In the output numpy arrays, which k-mer corresponds to a particular frequency in the array is not indicated. Still, k-mer order in the array is always the same for each time FastaToKmerFreq.py is run and for each sequence.  
To get the k-mer order in the arrays, the flag ```--get_kmers``` can be set. The program will just print to stdout the k-mer order for the specified size, ignoring all other arguments and terminating immediately after.  

**Example:**  

```python3 FastaToKmerFreq.py -kmer_size 3 -fasta_folder .\genomes\ -output_folder example_output -p 3 --get_kmers```  

**Output:**  

```AAA,AAC,AAG,AAT,ACA,ACC,ACG,ACT,AGA,AGC,AGG,ATA,ATC,ATG,CAA,CAC,CAG,CCA,CCC,CCG,CGA,CGC,CTA,CTC,GAA,GAC,GCA,GCC,GGA,GTA,TAA,TCA```  

The frequency value at index 0 of the output arrays corresponds to the frequency of "AAA" in the given sequence.
