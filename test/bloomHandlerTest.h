#include "gtest/gtest.h"
#include "bloom_filter.hpp"
#include "bloomHandler.hpp"

TEST(bloomIO, read_write){

  bloom_parameters parameters;
  parameters.projected_element_count    = 10000;
  parameters.false_positive_probability = 0.001;
  parameters.random_seed                =     1;
  parameters.compute_optimal_parameters()      ;

  bloomContainer created_blooms;

  for(int i = 0; i < 10; i++){
    created_blooms.add(parameters, 0);
    for(int j = 0 ; j < 100; j++){
      int r = rand(); 
      created_blooms.data[i]->bf.insert(r);
    }
  }
  
  std::string test = "test-index.bin";
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

  std::cerr << "HERE" << std::endl;
  
  for(int i = 0; i < 10; i++){
    ASSERT_EQ(created_blooms.data[i]->bf,
	      loaded_blooms.data[i]->bf);
  }  
};

