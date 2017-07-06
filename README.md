# fakeFastqGen

### How to compile
```
git clone https://github.com/hiraksarkar/fakeFastqGen.git
mkdir build
cmake ..
make
```
### How to run
```
src/handlefasta fastqgen -h
  fastq generator
Usage: build/src/handlefasta fastqgen [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -f,--fasta TEXT             the full path to the fasta file
  -o,--outdir TEXT            output directory
  -p,--paireend               library type
  -r,--read-len UINT          read length
```
### Example
```
src/handlefasta fastqgen -f 3_genomes.fasta -p -o test -r 20
```
create fake paired-end fastq files utilizing all of sequences present in the fasta file

### TODO
Have to add option to ignore short sequences
