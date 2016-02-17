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


extern "C"{
#include <zlib.h>
#include <stdint.h>
}

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)  

struct options{
   std::string file;
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
	globalOpts.file = optarg;
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
 Function input  :

 Function does   :

 Function returns:

*/
void sub()
{
}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
  int parse = parseOpts(argc, argv);
  
  if(globalOpts.file.empty()){
    std::cerr << "FATAL: could not open file: " << globalOpts.file << std::endl;
  }

  gzFile fp;
  kseq_t *seq;
  int l;

  fp = gzopen(globalOpts.file.c_str(), "r");
  seq = kseq_init(fp);

  while ((l = kseq_read(seq)) >= 0) {


    std::cerr << "len:" << seq->seq.l << std::endl;
    
    for(uint32_t i = 0; i < seq->seq.l -1; i++){

      uint64_t fkmer = 0;
      uint64_t rkmer = 0;
      int fail = 0;
      
      fail += dnaTobit(seq->seq.s, i, &fkmer);
      if(fail != 0){
	revcomp(fkmer, &rkmer);
	std::cerr << "fKMER:" << fkmer << std::endl;

	std::cerr << " expected: ";
	
	for(int j = i ; j < i + 32; j++){
	  std::cerr << seq->seq.s[j] ;
	}
	std::cerr << std::endl;
	std::cerr << " forward : ";
	printKmer(fkmer);
	std::cerr << " reverse : ";
	printKmer(rkmer);


      }
    }
  } 
  kseq_destroy(seq); 
  gzclose(fp); 
  
  return 0;
}
