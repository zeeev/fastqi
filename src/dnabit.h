#ifdef __cplusplus
extern "C" {
#endif

#ifndef DNABIT_H
#define DNABIT_H

#include <stdint.h>
#include <ctype.h>
#include <cstdio>
#include <inttypes.h>

  
#define BA 0
#define BT 1
#define BC 2
#define BG 3

  /* converts DNA to bits and packs it into a uint64 */
  
  bool dnaTobit( char * s, uint32_t start, uint64_t * results){
  
    uint32_t end          = start + 32 ;
    uint32_t stopShifting = end - 1    ; 
    
    for(uint32_t i = start; i < end; i++){
      uint8_t baseNum = toupper(s[i]);
      
      if(baseNum == 65){
	*results |= BA;
      }
      else if(baseNum == 84){
	*results |= BT;
      }
      else if(baseNum == 71){
	*results |= BG;
      }
      else if(baseNum == 67){
	*results |= BC;
      }
      else{
	return false;
      }
      
      if(i < stopShifting){
	*results <<= 2;
      }
    }
    return true;
  }

  uint8_t flip(uint8_t base){
    if(base == BG){
      return BC;
    }
    if(base == BC){
      return BG;
    }
    if(base == BA){
      return BT;
    }
    if(base == BT){
      return BA;
    }
    // not ideal
    return 0;
  }
  void revcomp(uint64_t seq, uint64_t * revcomp){
    uint64_t mask = 3;
    for(uint8_t i = 0; i <= 32; i++){
      *revcomp |= flip(seq & mask);
      if(i <= 30){
	seq >>= 2;
	*revcomp <<= 2;
      }
    }
  }

  void printKmer(uint64_t seq){

    uint64_t mask = 3;

    for(int i = 62; i >= 0; i-=2){

      if(((seq >> i) & mask) == BT){
	printf("T");	
      }
      else if(((seq >> i) & mask) == BA){
	printf("A");	
      }
      else if(((seq >> i) & mask) == BG){
	printf("G");	
      }
      else{
	printf("C");	
      }
    }

    
    printf("\n");
  }

  void printBin(uint64_t seq){

    uint64_t mask = 1;
    
    uint8_t res[64];
    
    for(int i = 63; i >= 0; i--){
      res[i] = (mask & seq);
      seq >>= 1;
    }

    for(int i = 0; i < 64 ; i++){
      printf("%u", res[i]);
    }
    printf("\n");
  }
  
  
#endif
#ifdef __cplusplus
}
#endif
