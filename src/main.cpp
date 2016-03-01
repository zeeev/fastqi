/*

This program was created at:  Wed Feb 17 08:20:38 2016
This program was created by:  Zev N. Kronenberg


Contact: zev.kronenber@gmail.com

Organization: University of Washington
    Seattle, Washington


The MIT License (MIT)

Copyright (c) <2016> <Zev N. Kronenberg>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


*/

#include <string>
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <map>
#include <vector>
#include <omp.h>

extern "C"{
#include <stdint.h>
}

#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "bloom_filter.hpp"
#include "bloomHandler.hpp"
#include "split.h"
#include "dnabit.h"

#define BLOCK_SIZE 1000

struct options{
  int nthreads;
  std::string file ;
  std::string index;
  std::string seqs ;
}globalOpts;

struct fastq_entry{
  kstring_t seqid;
  kstring_t seq  ;
  kstring_t sep  ;
  kstring_t qual ;
};


static const char *optString = "vhf:s:x:";

omp_lock_t lock;

void printHelp(void){
  std::cerr << std::endl;
  std::cerr << "Synopsis:                                          "<< std::endl;
  std::cerr << " Retrieve reads from fastq by kmer match (32-mer). "<< std::endl;
  std::cerr << std::endl;
  std::cerr << "Usage:                                             "<< std::endl;
  std::cerr << "        bgzip test.fq                              "<< std::endl;
  std::cerr << "        fqi -f test.fq.gz -s A,B,C,D > hits.fq     "<< std::endl;
  std::cerr << std::endl;
  std::cerr << "Required:  " << std::endl;
  std::cerr << " -f, <STRING> - A bgzipped fastq file.             "<< std::endl;
  std::cerr << " -s, <STRING> - A comma separated list of kmers.   "<< std::endl;
  std::cerr << "                The kmers must be 32 characters.   "<< std::endl;
  std::cerr << " -x, <INT>    - Number of threads [1]              "<< std::endl;
  std::cerr << "Output:  " << std::endl;
  std::cerr << " Fastq entries are printed to STDOUT               "<< std::endl;
}
//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
{
  int opt = 0;
  opt = getopt(argc, argv, optString);
  while(opt != -1){
    switch(opt){
    case 'h':
      {
	printHelp();
	break;
      }
    case 'x':
      {
	globalOpts.nthreads = atoi(((std::string)optarg).c_str());
	break;
      }
    case 'f':
      {
	globalOpts.file  = optarg;
	globalOpts.index = globalOpts.file + ".fqi";
	break;
      }
    case 's': 
      {
        globalOpts.seqs = optarg;
        break;
      }

    case '?':
      {
	break;
      }
    }    
    opt = getopt( argc, argv, optString ); 
  }
  return 1;
}
//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : a pointer to kstring_t

 Function does   : init 0

 Function returns: NA

*/

inline void initKstring(kstring_t * k){
  k->m = 0;
  k->l = 0;
  k->s = 0;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of kmers to search, vector of offsets,
                   bloomContainer

 Function does   : Goes through bloom container and finds offsets 
                   that need to be searched.

 Function returns: int ; > 0 good < 0 bad

*/

bool flatSearch(std::vector<uint64_t> & qKmers   ,
		std::map<uint64_t, bool> & of ,
		bloomContainer & bc            ){

  for(std::vector<uint64_t>::iterator it = qKmers.begin();
	it != qKmers.end(); it++){

    for(std::vector<bloomWrapper *>::iterator iz = bc.data.begin();
	iz != bc.data.end(); iz++){

      if((*iz)->bf.contains(*it)){
	of[(*iz)->fastqOffset] = true ;
      }
    }
  }

  return true;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : an offset position, vector<uint64_t> kmers to match, 
                   a string stream to load matches into

 Function does   : bool

 Function returns: int ; > 0 good < 0 bad

*/

bool getRecords(uint64_t pos, std::vector<uint64_t> & qKmers){

  BGZF * fp;

  fp = bgzf_open(globalOpts.file.c_str(), "r");

  if(bgzf_seek(fp, pos, 0) < 0){
    std::cerr << "FATAL: Could not seek within file." << std::endl;
    exit(1);
  }

  fastq_entry fq;
  initKstring(&fq.seqid);
  initKstring(&fq.seq  );
  initKstring(&fq.sep  );
  initKstring(&fq.qual );

  int nreads = 0;
  int state  = 1;

  uint64_t kmer = 0;
  
  while(state > 0 && nreads < BLOCK_SIZE){

    state = bgzf_getline(fp, '\n', &fq.seqid);
    state = bgzf_getline(fp, '\n', &fq.seq  );
    state = bgzf_getline(fp, '\n', &fq.sep  );
    state = bgzf_getline(fp, '\n', &fq.qual );    
    
    if(state < 0){
      break;
    }

    nreads++;

    std::map<uint64_t, bool> kmersInRead;
    
    for(uint32_t i = 0; i < fq.seq.l - 32 ; i++){
      kmer = 0;
      if(dnaTobit(fq.seq.s, i, &kmer)){
	kmersInRead[kmer] = true;
      }
    }

    bool hit = false;
    
    for(std::vector<uint64_t>::iterator iz = qKmers.begin();
	iz != qKmers.end(); iz++){

      if(kmersInRead.find(*iz) != kmersInRead.end()){
	hit = true;
      }
    }

    if(hit){
      omp_set_lock(&lock);

      std::cout << fq.seqid.s << "\n" << fq.seq.s
	 << "\n" << fq.sep.s << "\n" << fq.qual.s << "\n";
      
      omp_unset_lock(&lock);
    }
  }

  bgzf_close(fp);
  return true;

}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : an offset position

 Function does   : processes a chunk
                   
 Function returns: int ; > 0 good < 0 bad

*/
int processChunk(uint64_t pos, bloomWrapper * bfw){

  BGZF * fp;

  std::vector<uint64_t> kmers;
  
  fp = bgzf_open(globalOpts.file.c_str(), "r");

  if(bgzf_seek(fp, pos, 0) < 0){
    std::cerr << "FATAL: Could not seek within file." << std::endl;
    exit(1);
  }

  int nreads = 0;
  int state  = 1;

  uint64_t kmer = 0;
  fastq_entry fq;

  initKstring(&fq.seqid);
  initKstring(&fq.seq  );
  initKstring(&fq.sep  );
  initKstring(&fq.qual );

  while(state > 0 && nreads < BLOCK_SIZE){

    state = bgzf_getline(fp, '\n', &fq.seqid);
    state = bgzf_getline(fp, '\n', &fq.seq  );
    state = bgzf_getline(fp, '\n', &fq.sep  );
    state = bgzf_getline(fp, '\n', &fq.qual );

    if(state < 0){
      break;
    }

    for(uint32_t i = 0; i < fq.seq.l - 32 ; i++){
      kmer = 0;
      if(dnaTobit(fq.seq.s, i, &kmer)){
	kmers.push_back(kmer); 
	
      }
    }
    nreads += 1;
  }

  for(int i = 0; i < kmers.size(); i++){
    bfw->bf.insert(kmers[i]);
  }
  
  bgzf_close(fp);

  return 1;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : a vector for compressed index positions (tell)

 Function does   : loads up the offsets for each block (100 reads) 
                   starting at read

 Function returns: int ; > 0 good < 0 bad

*/
int getReadOffsets(std::vector<uint64_t> & offsets)
{

  BGZF * fp;
    
  fp = bgzf_open(globalOpts.file.c_str(), "r");
  
//  if(bgzf_index_load(fp, globalOpts.file.c_str(), ".gzi") <0){
//    std::cerr << "FATAL: Could not open index file.";
//    exit(1);
//  }
  
  uint64_t nreads = 0;
  fastq_entry fq;

  initKstring(&fq.seqid);
  initKstring(&fq.seq  );
  initKstring(&fq.sep  );
  initKstring(&fq.qual );
  
  int state = 1;
  
  while (state > 0 ){
    uint64_t tellp = bgzf_tell(fp)           ;

    state = bgzf_getline(fp, '\n', &fq.seqid);
    state = bgzf_getline(fp, '\n', &fq.seq  );
    state = bgzf_getline(fp, '\n', &fq.sep  );
    state = bgzf_getline(fp, '\n', &fq.qual );

    if( (nreads % BLOCK_SIZE) == 0 || nreads == 0){
      offsets.push_back(tellp);
    }
    nreads+= 1;
  }
  
  if(bgzf_close(fp) < 0){
    std::cerr << "FATAL: problem closing file handle" << std::endl;
    exit(1);
  }

  
  if(offsets.size() > 0){
    return 1;
  }
  else{
    return 0;
  }  
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector to load

 Function does   : turns a string into a uint64_t;

 Function returns: bool

*/


bool parseSeqs(std::vector<uint64_t> & toLoad){
  

  std::vector<std::string> seqs = split(globalOpts.seqs, ",");
  
  for(std::vector<std::string>::iterator it = seqs.begin();
      it != seqs.end(); it++){

  
    if((*it).size() != 32){
      std::cerr << "WARNING: skipping kmer that is not 32 bases." << std::endl;
      continue;
    }
 
    uint64_t kmer = 0;
    
    char * chr = strdup((*it).c_str());

    if(dnaTobit(chr, 0, &kmer)){
      toLoad.push_back(kmer);
    }
  }
    std::cerr << "INFO: going to search for " 
	      << toLoad.size() << " kmers." << std::endl;

  return true;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : nothing

 Function does   : the ugly parts of building the index

 Function returns: bool

*/

bool buildIndex(void){

  bloomIO controller(globalOpts.index);
  if(controller.openForRead()){
    return true;
  }

  std::vector<uint64_t> offsets;
  if( getReadOffsets(offsets) < 0 ){
    std::cerr << "FATAL: problem reading offsets" << std::endl;
    exit(1);
  }
  
  bloom_parameters parameters;
  parameters.projected_element_count    = 5000  ;
  parameters.false_positive_probability = 0.001 ;
  parameters.random_seed                =rand() ;
  parameters.compute_optimal_parameters()       ;
  
  bloomContainer created_blooms;
  
  int i = 0;
  for(i = 0; i < offsets.size(); i++){
    created_blooms.add(parameters, offsets[i]);
  }
  
#pragma omp parallel for schedule(dynamic, 3)
  for(i = 0; i < offsets.size(); i++){
    processChunk(offsets[i], created_blooms.data[i]);

    if((i % 1000) == 0){
      std::cerr << "INFO: processed " << i
		<< " chunks of " << offsets.size() << std::endl;
    }  
  }

  bloomIO controller2(globalOpts.index);
  
  controller.openForWrite();
    
  controller.write(created_blooms);
 
  return true;
}
//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
  
  globalOpts.nthreads = 1;

  int parse = parseOpts(argc, argv);

  omp_set_num_threads(globalOpts.nthreads);

  if(globalOpts.file.empty()){
    std::cerr << "FATAL: no file (-f) " << std::endl;
    printHelp();
    exit(1);
  }

  if(globalOpts.seqs.empty()){
    std::cerr << "FATAL: no seqs (-s) " << std::endl;
    printHelp();
    exit(1);
  }

  std::vector<uint64_t> toFind ;

  if(!parseSeqs(toFind)){
    std::cerr << "FATAL: problem loading query sequences." << std::endl;
  }

  buildIndex();

  bloomContainer loaded_blooms;
  
  {
    bloomIO controller(globalOpts.index);
    controller.openForRead();    
    controller.read(loaded_blooms);
  }

  std::map<uint64_t, bool> offsetsToRead;
    
  flatSearch(toFind,offsetsToRead,loaded_blooms);

  std::cerr << "There are " << offsetsToRead.size()
	    << " file blocks that must be read out of " 
	    << loaded_blooms.data.size() << std::endl;

  std::map<uint64_t, bool>::iterator it = offsetsToRead.begin();
  
#pragma omp parallel for schedule(dynamic, 3)    
  for(uint32_t i = 0; i < offsetsToRead.size(); i++){
    getRecords(it->first, toFind);
    it++;
  }
    
  return 0;
}
