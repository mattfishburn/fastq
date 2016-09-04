# fastq
kmer counting of FASTQ files

# Notes

I did not test any unicode compatibility of the file.  Assumed ASCII.

If there are ties in the output count, the order things are output is random.  So, if you are printing the top 30, and there is a five-way tie at count 100 for 29th place, the values output for 29th and 30th are randomly two of the five.

I treated the letter 'N' as though it were it's own letter.  I would be interested in what a biologist would want to do with that letter.

Took about three hours:
* About an hour to get the python implementation working and the c++ implementation designed based on
what I learned from python.
* Another house to code up the c++ implementation (in Visual Studio).
* Another hour to get the code to compile on ubuntu

I profiled the VS C++ implementation; about half of the implementation
time is spent checking if the kmer already exists in the hashmap
countFromKmer along with inserting it if it does not.  The next step I
would take is to modify the value allocator for unordered_map to
return zero, and just use ++countFromKmer[sub] instead of checking if
the key exists and inserting it if it doesn't.  This would save a hash
and a lookup.  After that, I would probably multithread the problem
using a dataflow model with a single producer consumer.  Another
performance enhancement would be to improve the substring pulls.

# Approach 

## python reference implementation

Started with a python implementation in the python folder, python
implementation took about 7.5 seconds to run (Windows), breakdown was
100 ms reading and parsing file, 5000 ms getting kmers from file lines
and histogramming, and then 2200 milliseconds sorting histogram
contents.

Most of the kmers are unique in the file I tested.  There are also 'N' letters (!), not listed on wikipedia page, other sources imply these letters are when that DNA letter could not be read.

## Visual Studio c++ implementation

Uses a hashmap / dictionary to map from kmer to count.  Minimized sort
time using a priority queue to track only the topcount biggest values
when outputting the most seen kmers.

Would expect the runtime of the various parts to be:

* File read: O(size) where size is the filesize
* Histogramming: O(uniques) where uniques is the total number of kmers (using a hashmap / dictionary)
* Topcount: O(uniques + topcount * (lg topcount) * (lg uniques)) where topcount is the number to output (probabilistic)

## c++ implementation (linux)

Clang wouldn't compile the decltype statement with the tuple, I
changed it to use a standard comparator.