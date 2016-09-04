# fastq
kmer counting of FASTQ files

Compiles on an ubuntu machine with clang-3.5 installed.  No external libraries are required.  To run:
```
make all
./kmercount fastq-filename kmersize topcount
```

There is an example file in the example folder, so you should be able to run:
```
./kmercount example/ERR055763.filt.fastq 30 25
```

# Implementation Overview

The program has three main parts:
* Read input from file into vector lines
* Histogram k-mers from the lines vector using a string->int hashmap (unordered_map)
* After histograming, use a priority queue to track the top kmers, and output those top k-mers and the count to the command line.

# Notes

I did not test any unicode compatibility of the file.  Assumed ASCII.

If there are ties in the output count, the order kmers are output is random.  So, if you are printing the top 30, and there is a five-way tie at count 100 for 29th place, the values output for 29th and 30th are randomly two of the five.

Program exits with error code 1 when usage is incorrect, or when fewer than topcount kmers are found.

I treated the letter 'N' in the fastq file as though it were it's own letter.  I would be interested in how a biologist would  interpret the letter.

Overall, took about three hours:
* About an hour to get the python implementation working and the c++ implementation designed based on what I learned from python.
* Another house to code up the c++ implementation (in Visual Studio).
* Another hour to get the code to get ubuntu setup with the right packages and the code to compile.

I tested the VS c++ and ubuntu c++ against python on one of the files.

I profiled the VS C++ implementation; about a quarter of the
time is spent checking if the kmer already exists in the hashmap
countFromKmer, and about another quarter is  inserting it if it does not.  The next optimization I
would take is to modify the value allocator for unordered_map to
return zero, and just use ++countFromKmer[sub] instead of checking if
the key exists along with inserting it if it doesn't (lines 80 - 86).  This would save a hash
and a lookup.  After that, I would probably multithread the problem
using a dataflow model with a single producer consumer.  Other
performance enhancement would be to improve the substring pulls, improving the hash speed, trying a better hash_map, and more multi-threading like a mapReduce algorithm.

# Approach 

## python reference implementation

Started with a python implementation in the python folder, python
implementation took about 7.5 seconds to run (Windows), breakdown was
100 ms reading and parsing file, 5000 ms getting kmers from file lines
and histogramming, and then 2200 milliseconds sorting histogram
contents.

Most of the kmers are unique in the file I tested.  There are also 'N' letters (!), not listed on wikipedia page, other sources imply these letters signify that DNA letter could not be read.  Handled these as letter 'N'.

## Visual Studio c++ implementation

Uses a hashmap / dictionary to map from kmer to count.  Minimized sort
time using a priority queue to track only the topcount biggest values
when outputting the most seen kmers.  Takes a little over half the time that python takes.

Would expect the runtime of the various parts to be:

* File read: O(size) where size is the filesize
* Histogramming: O(uniques) where uniques is the total number of kmers (using a hashmap / dictionary)
* Topcount: O(uniques + topcount * (lg topcount) * (lg uniques)) where topcount is the number to output (probabilistic)

## c++ implementation (linux)

Clang wouldn't compile the decltype statement with the tuple, I
changed it to use a standard comparator.
