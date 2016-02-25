#ifndef BLOOMHANDLER_H
#define BLOOMHANDLER_H

#include <fstream>
#include <stdint.h>
#include "bloom_filter.hpp"

class bloomWrapper{
  friend class bloom_filter     ;
  friend class bloom_parameters ;
  
public:

  bloom_filter                bf ;
  uint32_t                    id ;
  uint32_t         parentPointer ;
  uint64_t           fastqOffset ;
  std::vector<uint32_t> children ;

  bloomWrapper(const bloom_parameters& p, uint64_t of) :
    bf(p), parentPointer(0), fastqOffset(of)
  {} 
};


class bloomContainer{
public:
  
  friend class     bloom_filter ;
  friend class bloom_parameters ;
  friend class     bloomWrapper ;
  
  std::vector<bloomWrapper *> data;

  ~bloomContainer()
  {
    for(std::vector<bloomWrapper *>::iterator it = data.begin();
	it != data.end(); it++){
      delete (*it);
    }
  }

  inline void add(const bloom_parameters& p, uint64_t of){
    bloomWrapper * bf = new bloomWrapper(p, of);
    data.push_back(bf);
  }
  
};

class bloomIO{
  friend class     bloom_filter ;
  friend class bloom_parameters ;
  friend class     bloomWrapper ;

private:
  uint32_t          lucky ;
  bool         write_lock ;
  std::string    fileName ;
  std::fstream         fs ;

public:

  bloomIO(std::string & fn) : fileName(fn), write_lock(false)
  {
    lucky = 1984;
  }

  ~bloomIO()
  {
    fs.close();
  }
  
  inline bool openForRead(void){
        fs.open(fileName.c_str(), std::fstream::in | std::fstream::binary);

	uint64_t magic;
	
	fs.read((char*)&magic, sizeof(uint64_t));

	if(magic == lucky){
	  std::cerr << "INFO: good to read index" << std::endl;
	}
	
	return true;
  }

  inline bool openForWrite(void){

    if(fs.is_open()) return true;

    fs.open(fileName.c_str(), std::fstream::out | std::fstream::binary);

    if(fs.is_open()){
      fs.write((char *) &lucky, sizeof(uint64_t));
      return true;
    }
    else{
      return false;
    }
  }

  inline bool write(std::vector<bloomWrapper *> bws){

    if(write_lock) return false;

    if(!fs.is_open()){
      openForWrite();
    }

    for(std::vector<bloomWrapper *>::iterator it = bws.begin();
	it != bws.end(); it++){
      fs.write((char *)(*it), sizeof(bloomWrapper));
    }

    write_lock = true;
    return true;
  }  
};



#endif
