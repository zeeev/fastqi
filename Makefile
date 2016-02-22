######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CXX=g++

fqi: bin libhts.a
	cd bin && $(CXX) -L ../htslib -I ../htslib -I ../src/ -lz ../src/main.cpp ../htslib/libhts.a  -o fqi

bin:
	mkdir bin

libhts.a:
	cd htslib && make

