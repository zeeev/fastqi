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
#include "dnabit.h"
#include <map>
#include <vector>

extern "C"{
#include <stdint.h>
}

#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "bloom_filter.hpp"
#include "bloomHandler.hpp"

struct options{
  std::string file ;
  std::string bft  ;
  std::string index;
}globalOpts;

struct fastq_entry{
  kstring_t seqid;
  kstring_t seq  ;
  kstring_t sep  ;
  kstring_t qual ;
};


static const char *optString = "hf:";

//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
{
  int opt = 0;
  opt = getopt(argc, argv, optString);
  while(opt != -1){
    switch(opt){
    case 'h':
      {
	break;
      }
    case 'f':
      {
	globalOpts.file  = optarg;
	globalOpts.bft   = globalOpts.file + ".bft";
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
  
  while(state > 0 && nreads < 100){

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
      std::cout << fq.seqid.s << "\n" << fq.seq.s
	 << "\n" << fq.sep.s << "\n" << fq.qual.s << "\n";
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

  std::cerr << "bin offset: " << pos << std::endl;
  
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

  while(state > 0 && nreads < 100){

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

  std::cerr << "N kmers in bin: " << kmers.size() << std::endl;

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

    if( (nreads % 100) == 0 || nreads == 0){
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

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
  int parse = parseOpts(argc, argv);
  
  if(globalOpts.file.empty()){
    std::cerr << "FATAL: no file (-f) " << std::endl;
    exit(1);
  }

  std::vector<uint64_t> offsets;
  if( getReadOffsets(offsets) < 0 ){
    std::cerr << "FATAL: problem reading offsets" << std::endl;
    exit(1);
  }

  bloom_parameters parameters;
  parameters.projected_element_count    = 10000;
  parameters.false_positive_probability = 0.001; 
  parameters.random_seed                =     1;
  parameters.compute_optimal_parameters()      ;

  bloomContainer created_blooms;

  int i = 0;
  for(i = 0; i < offsets.size(); i++){
    created_blooms.add(parameters, offsets[i]);
  }
  
  for(i = 0; i < offsets.size(); i++){
    processChunk(offsets[i], created_blooms.data[i]);
  }

  std::string test = "testi.bin";
  {
    bloomIO controller(test);
    
    controller.openForWrite();

    controller.write(created_blooms);
  }

  bloomContainer loaded_blooms;
  
  {
    bloomIO controller(test);
    
    controller.openForRead();
    
    controller.read(loaded_blooms);
  }

  std::vector<uint64_t> toFind ;
  std::map<uint64_t, bool> offsetsToRead;
  
  uint64_t k1 = 0;
  uint64_t k2 = 14077770871635063417;
  uint64_t k3 = 14077770871635063417;
  uint64_t k4 = 14077770871635063417;
  
  toFind.push_back(k1);
  toFind.push_back(k2);
  toFind.push_back(k3);
  toFind.push_back(k4);
  
  flatSearch(toFind,offsetsToRead,loaded_blooms);

  std::cerr << "There are " << offsetsToRead.size()
	    << " file offsets that must be read out of " 
	    << loaded_blooms.data.size() << std::endl;
    
  
  for(std::map<uint64_t, bool>::iterator it = offsetsToRead.begin();
      it != offsetsToRead.end(); it++){
    getRecords(it->first, toFind);
  }
    
  return 0;
}
