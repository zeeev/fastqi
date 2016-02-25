#ifndef BLOOMHANDLER_H
#define BLOOMHANDLER_H

#include <fstream>
#include <stdint.h>
#include "bloom_filter.hpp"

class bloomWrapper{
  friend class bloom_filter     ;
  
public:

  bloom_filter                bf ;
  uint32_t                    id ;
  uint32_t          parentArrayI ;
  uint64_t           fastqOffset ;
  std::vector<uint32_t> children ;

  bloomWrapper(const bloom_parameters& p, uint64_t of) :
    bf(p), parentArrayI(0), fastqOffset(of), id(0)
  {}

  inline bool write(std::fstream & fs){
    if(!fs.is_open()) return false;

    // writing number of salts
    fs.write((char *) &bf.salt_count_, sizeof(unsigned int) );

    // writing salt vector
    for(unsigned int i = 0; i < bf.salt_count_ ; i++){
      fs.write((char *) &bf.salt_[i], sizeof(bloom_type) );
    }

    // writing table size
    fs.write((char *) &bf.table_size_,
	     sizeof(unsigned long long int) );

    // writing raw_table_size
    fs.write((char *) &bf.raw_table_size_,
	     sizeof(unsigned long long int) );
    
    // writing table
    for(unsigned long long int i = 0;
	i < bf.salt_count_ ; i++){
      fs.write((char *) &bf.bit_table_[i],
	       sizeof(unsigned char) );
    }

    // writing projected element count
    fs.write((char *) &bf.projected_element_count_,
	     sizeof(unsigned long long int) );

    // writing inserted_element_count_
    fs.write((char *) &bf.projected_element_count_,
	     sizeof( unsigned int ) );

    // writing random seed
    fs.write((char *) &bf.random_seed_,
	     sizeof( unsigned long long int ) );

    // prob
    // writing random seed
    fs.write((char *) &bf.desired_false_positive_probability_,
	     sizeof( double ) );

    return true;
  }

  
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
};



#endif
