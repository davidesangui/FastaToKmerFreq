# FastaToKmerFreq
Simple script to reduce DNA sequences (in Fasta format) to K-mer frequencies.
## Installation
To use FastaToKmerFreq, simply download the file FastaToKmerFreq.py from this repository and run it with Python3.
The script requires numpy, which can be obtained via pip:
```
pip install numpy
```
## How it works
FastaToKmerFreq takes as input a folder of fasta files. All the contigs in each fasta file are reduced to a single set of k-mer frequencies representing the entire sequence stored in the the fasta file; alternatively, each contig is reduced to a contig-specific set of k-mer frequencies. The k-mer size is user-defined.  
In theory, there are 4^kmer_size possible kmers. But when considering forward and reverse k-mers and palindromic k-mers, there is a lower number of distinct k-mers (e.g with k=4, there are 136 distinct kmers; for more details: https://stackoverflow.com/questions/40952719/algorithm-to-collapse-forward-and-reverse-complement-of-a-dna-sequence-in-python).  
FastaToKmerFreq reduces a sequence to a set of distinct k-mer frequencies by counting the k-mers in both the forward and reverse sequence.  
All k-mers containing characters other than {'A','C','G','T'} are ignored.  
K-mer frequencies are saved as numpy arrays. Depending on the output type requested, the output is one numpy array per fasta file, one numpy matrix with as many rows as fasta files, or one numpy matrix with as many rows as contigs in all fasta files.  
When multiple cores are available, different Fasta files can be reduced in parallel.
## Usage
FastaToKmerFreq.py requires three mandatory arguments:
- ```-kmer_size``` Size of k-mers. Since all possible k-mers are stored in memory at once, increasing the k-mer size quickly increases memory consumption. A rough estimation of memory consumption (in bytes) can be performed with the following formula: kmer_size * numCores * (4^kmer_size)/2
- ```-fasta_folder``` Path to folder containing (only) the fasta files to reduce to kmer frequencies.
- ```-output_folder``` Path to output folder to store kmer frequencies and index file. Created if missing.

**Example:**  

```python3 FastaToKmerFreq.py -kmer_size 3 -fasta_folder .\genomes\ -output_folder example_output -p 3``` 

```-p``` defines the number of processes for parallelization. Setting a number greater than the number of fasta files is useless.  

### Output type  
```-output_type``` option allows to define the output type desired. Three different types are available:
- single_files: (default) Each fasta file is entirely reduced to a single set of k-mer frequencies. The output consists of one numpy array per fasta file.  
In this repository, the examples are: example_output/Ecoli.fna.npy, example_output/Paeruginosa.fna.npy, example_output/Senterica.fna.npy.
- fasta_matrix: Each fasta file is entirely reduced to a single set of k-mer frequencies. The output consists of a single numpy matrix where each row stores the k-mer frequency of the corresponding fasta file. The fasta_matrix_index.txt file lists all fasta files in the same order of appearance in the matrix. The fasta file listed in the first line of the index file corresponds to the frequency set at row 0 in the output matrix, and so on.  
In this repository, the examples are: example_output/fasta_matrix.npy, example_output/fasta_matrix_index.txt
- contig_matrix: Each contig in all fasta files is reduced to a set of k-mer frequencies. The output consists of a single numpy matrix where each row stores the k-mer frequency of the corresponding contig. The contig_matrix_index.txt file lists all contigs in the same order of appearance in the matrix. The contig listed in the first line of the index file corresponds to the frequency set at row 0 in the output matrix, and so on.  
In this repository, the examples are: example_output/contig_matrix.npy, example_output/contig_matrix_index.txt

**Example:**  

```python3 FastaToKmerFreq.py -kmer_size 3 -fasta_folder .\genomes\ -output_folder example_output -p 3 -output_type contig_matrix``` 

### K-mer order
In the output numpy arrays, which k-mer corresponds to a particular frequency in the array is not indicated. Nevertheless, k-mer order in the array is always the same for each time FastaToKmerFreq.py is run and for each sequence.  
To get the k-mer order in the arrays, the flag ```--get_kmers``` can be set. The program will just print to stdout the k-mer order for the specified size, ignoring all other arguments and terminating immediately after. 

**Example:**  

```python3 FastaToKmerFreq.py -kmer_size 3 -fasta_folder .\genomes\ -output_folder example_output -p 3 --get_kmers```  

**Output:**  

```AAA,AAC,AAG,AAT,ACA,ACC,ACG,ACT,AGA,AGC,AGG,ATA,ATC,ATG,CAA,CAC,CAG,CCA,CCC,CCG,CGA,CGC,CTA,CTC,GAA,GAC,GCA,GCC,GGA,GTA,TAA,TCA```  

The frequency value at column 0 of the output arrays corresponds to the frequency of "AAA" in the given sequence.
