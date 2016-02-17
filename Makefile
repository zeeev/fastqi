######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CXX=g++

fqi: bin
	cd bin && $(CXX) -I ../src/ -lz ../src/main.cpp -o fqi

bin:
	mkdir bin
