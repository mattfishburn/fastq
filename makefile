all: kmercount

kmercount: ubuntu/fastq.cpp
	cd ubuntu && (clang++-3.5 fastq.cpp --std=c++11)
	cp ubuntu/a.out kmercount
