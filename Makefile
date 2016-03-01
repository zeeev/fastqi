######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CXX=c++ -lstdc++

fqi: bin libhts.a runtest
	cd bin && $(CXX) -g -O0 -L ../htslib -I ../bloom -I ../htslib -I ../src/ -lz ../src/main.cpp ../src/split.cpp ../htslib/libhts.a  -lpthread -o fqi

bin:
	mkdir bin
libhts.a:
	cd htslib && make
runtest: test/mainTest
	./test/mainTest
test/mainTest: gtest-1.7.0/build/libgtest.a
	$(CXX) -I bloom -I src -I htslib -I gtest-1.7.0/include -L gtest-1.7.0/build -lpthread -lgtest test/main.cpp gtest-1.7.0/build/*.a -o test/mainTest

gtest-1.7.0/build/libgtest.a:
	cd gtest-1.7.0 && mkdir build && cd build && cmake .. && make
clean:
	rm test*bin && rm test/mainTest && rm test-index.bin
