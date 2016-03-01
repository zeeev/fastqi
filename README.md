# FASTQI 
### The fastqi repo

The repo will contain tools for indexing fastq files by kmer.
Currently the main tool is under development and not ready to use.

However, the google bloom filter wrapper is ready.

### bloomHandler.hpp

Contains classes to wrap googles bloom filter

 - [X] Serialize vector of bloom_filters
 - [ ] build bloom filter tree

### fqi

Fqi builds a fastq.gz index based on bloom filters.
Currently it is a flat index divided by file chunk.
Soon fqi will have a bloom tree for faster retrieval.

#### Building fqi
```
git clone --recursive https://github.com/zeeev/fastqi.git
cd fastqi
make
```

#### Running fqi

```
Synopsis:
 Retrieve reads from fastq by kmer match (32-mer).

Usage:
        bgzip test.fq
        fqi -f test.fq.gz -s A,B,C,D > hits.fq

Required:
 -f, <STRING> - A bgzipped fastq file.
 -s, <STRING> - A comma separated list of kmers.
                The kmers must be 32 characters.
Output:
 Fastq entries are printed to STDOUT
```







