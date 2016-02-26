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
  bloomWrapper(void) :
    parentArrayI(0), fastqOffset(0), id(0)
  {}

  inline bool read(std::fstream & fs){
    if(!fs.is_open()) return false;

    fs.read((char *) &fastqOffset,
	     sizeof(uint64_t) );

    fs.read((char *) &bf.salt_count_,
	     sizeof(unsigned int) );

    fs.read((char *) &bf.table_size_,
	     sizeof(unsigned long long int) );

    fs.read((char *) &bf.raw_table_size_,
	     sizeof(unsigned long long int) );

    fs.read((char *) &bf.projected_element_count_,
	     sizeof(unsigned long long int) );

    fs.read((char *) &bf.inserted_element_count_,
	     sizeof(unsigned int) );

    fs.read((char *) &bf.random_seed_,
	     sizeof(unsigned long long int) );

    fs.read((char *) &bf.desired_false_positive_probability_,
	     sizeof(double) );

    unsigned int v;
    
    for(unsigned int i = 0; i < bf.salt_count_ ; i++){
      fs.read((char *) &v, sizeof(unsigned int) );
      bf.salt_.push_back(v);
    }

    delete[] bf.bit_table_;
    bf.bit_table_ = new unsigned char[static_cast<std::size_t>(bf.raw_table_size_)];
    
    for(unsigned long long int i = 0;
	i < bf.table_size_ ; i++){
      fs.read((char *) &bf.bit_table_[i],
	       sizeof(unsigned char) );
    }
    return true;
  }
    
  inline bool write(std::fstream & fs){
    if(!fs.is_open()) return false;

    fs.write((char *) &fastqOffset,
	     sizeof(uint64_t) );
    
    fs.write((char *) &bf.salt_count_,
	     sizeof(unsigned int) );

    fs.write((char *) &bf.table_size_,
	     sizeof(unsigned long long int) );

    fs.write((char *) &bf.raw_table_size_,
	     sizeof(unsigned long long int) );

    fs.write((char *) &bf.projected_element_count_,
	     sizeof(unsigned long long int) );

    fs.write((char *) &bf.inserted_element_count_,
	     sizeof(unsigned int) );

    fs.write((char *) &bf.random_seed_,
	     sizeof(unsigned long long int) );

    fs.write((char *) &bf.desired_false_positive_probability_,
	     sizeof(double) );
    
    for(unsigned int i = 0; i < bf.salt_count_ ; i++){
      fs.write((char *) &bf.salt_[i], sizeof(unsigned int) );
    }

    for(unsigned long long int i = 0;
	i < bf.table_size_ ; i++){
      fs.write((char *) &bf.bit_table_[i],
	       sizeof(unsigned char) );
    }
    return true;
  }

  
};


class bloomContainer{
public:
    
  std::vector<bloomWrapper *> data;

  ~bloomContainer()
  {
    for(std::vector<bloomWrapper *>::iterator it = data.begin();
	it != data.end(); it++){
      delete (*it);
    }
  }

  inline void add(void){
     bloomWrapper * bf = new bloomWrapper;
     this->data.push_back(bf);
  }
  
  inline void add(const bloom_parameters& p, uint64_t of){
    bloomWrapper * bf = new bloomWrapper(p, of);
    this->data.push_back(bf);
  }
  
};

class bloomIO{

private:
  uint32_t          lucky ;
  bool         write_lock ;
  bool          read_lock ;
  std::string    fileName ;
  std::fstream         fs ;

public:

  bloomIO(std::string & fn) : fileName(fn), write_lock(false), read_lock(false)
  {
    lucky = 1984;
  }

  ~bloomIO()
  {
    fs.close();
  }

  inline bool read(bloomContainer & blooms){
    if(!fs.is_open()) return false;
    if(read_lock)     return false;

    uint64_t nFrozenBlooms;
    
    this->fs.read((char *) &nFrozenBlooms, sizeof(uint64_t));

    std::cerr << "Going to read: " << nFrozenBlooms << std::endl;
    
    return true;
  }

  
  inline bool write(bloomContainer & blooms){
    if(!fs.is_open()) return false;
    if(write_lock)    return false;
    uint64_t nblooms = blooms.data.size();
    
    this->fs.write((char *) &nblooms, sizeof(uint64_t));

    std::cerr << "Going to write: " << nblooms << std::endl;
    
    for(std::vector<bloomWrapper *>::iterator it = blooms.data.begin();
	it != blooms.data.end(); it++){
      (*it)->write(fs);
    }
    
    write_lock = true;
    return true;
  }
  
  inline bool openForRead(void){
        this->fs.open(fileName.c_str(), std::fstream::in
		      | std::fstream::binary);

	uint64_t magic;
	
	this->fs.read((char*)&magic, sizeof(uint64_t));

	if(magic == lucky){
	  std::cerr << "INFO: good to read index" << std::endl;
	}
	
	return true;
  }

  inline bool openForWrite(void){

    if(this->fs.is_open()) return true;

    this->fs.open(fileName.c_str(), std::fstream::out
		  | std::fstream::binary);

    if(fs.is_open()){
      this->fs.write((char *) &lucky, sizeof(uint64_t));
      return true;
    }
    else{
      return false;
    }
  }
};



#endif
