######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CXX=c++ -lstdc++

fqi: bin libhts.a runtest
	cd bin && $(CXX) -O3 -L ../htslib -I ../bloom -I ../htslib -I ../src/ -lz ../src/main.cpp ../htslib/libhts.a  -o fqi

bin:
	mkdir bin
libhts.a:
	cd htslib && make
runtest: test/mainTest
	./test/mainTest
test/mainTest: gtest-1.7.0/build/libgtest.a
	$(CXX) -I bloom -I src -I htslib -I gtest-1.7.0/include -L gtest-1.7.0/build -lgtest test/main.cpp -o test/mainTest

gtest-1.7.0/build/libgtest.a:
	cd gtest-1.7.0 && mkdir build && cd build && cmake .. && make
clean:
	rm *bin && rm test/mainTest && rm test-index.bin
