######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CXX=c++ -lstdc++

fqi: bin libhts.a
	cd bin && $(CXX) -O3 -L ../htslib -I ../bloom -I ../htslib -I ../src/ -lz ../src/main.cpp ../htslib/libhts.a  -o fqi

bin:
	mkdir bin

libhts.a:
	cd htslib && make

