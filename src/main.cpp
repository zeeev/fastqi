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



// This should overwrite the stupid opaque pointer

struct options{
  std::string file;
  std::string index;
}globalOpts;

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
 Function input  : an offset position

 Function does   : processes a chunk
                   

 Function returns: int ; > 0 good < 0 bad

*/
int processChunk(uint64_t pos){

  BGZF * fp;

  std::vector<uint64_t> kmers;

  std::cerr << "bin offset: " << pos << std::endl;
  
  fp = bgzf_open(globalOpts.file.c_str(), "r");

  if(bgzf_index_load(fp, globalOpts.file.c_str(), ".gzi") <0){
    std::cerr << "FATAL: Could not open index file.";
    exit(1);
  }

  if(bgzf_seek(fp, pos, 0) < 0){
    std::cerr << "FATAL: Could not seek within file." << std::endl;
    exit(1);
  }

  kstring_t fastq_line  ;  

  fastq_line.m = 0;
  fastq_line.l = 0;
  fastq_line.s = 0;
    
  uint8_t  fqRecord = 0;
  uint64_t nreads   = 0;
  uint64_t kmer     = 0;

 
  int state = 1;

  while (state > 0){

    if(nreads > 99){
      break;
    }
    
    state = bgzf_getline(fp, '\n', &fastq_line);

    if(state < 0){
      break;
    }
    fqRecord  +=1;
    if(fqRecord > 4){
      fqRecord = 1;
    }

    if(fqRecord == 1){
      nreads += 1;

      std::cerr << fastq_line.s << std::endl;

      
      for(uint32_t i = 0; i < fastq_line.l ; i++){
	kmer = 0;
	if(dnaTobit(fastq_line.s, i, &kmer)){
	  kmers.push_back(kmer);
	}
      }
    }
  }

  std::cerr << "N kmers in bin: " << kmers.size() << std::endl;
  
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
  
  if(bgzf_index_load(fp, globalOpts.file.c_str(), ".gzi") <0){
    std::cerr << "FATAL: Could not open index file.";
    exit(1);
  }
  
  kstring_t fastq_line ;

  fastq_line.m = 0;
  fastq_line.l = 0;
  fastq_line.s = 0;
  
  uint8_t   fqRecord = 0;
  uint64_t  nreads   = 0;

  int state = 1;
  
  while (state > 0){

    uint64_t tellp = bgzf_tell(fp)             ;
    state = bgzf_getline(fp, '\n', &fastq_line) ;

    if(state < 0){
      break;
    }
    fqRecord  +=1;
    if(fqRecord > 4){
      fqRecord = 1;
    }
    if(fqRecord == 2){
      nreads += 1;    
      if((nreads % 100) == 0 || nreads == 1){
	offsets.push_back(tellp);
      }
    }
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
  if( getReadOffsets(offsets) < 0 ){
    std::cerr << "FATAL: problem reading offsets" << std::endl;
    exit(1);
  }  
  for(std::vector<uint64_t>::iterator it = offsets.begin();
      it != offsets.end(); it++){
    processChunk(*it);
  }
  
  return 0;
}
