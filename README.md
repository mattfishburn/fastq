# fastq
kmer counting of FASTQ files

# Notes

I did not test any unicode compatibility of the file.  Assumed ASCII.

I treated the letter 'N' as though it were it's own letter.  I would be interested in what a biologist would want to do with that letter.

# python reference implementation
Started with a python implementation in the python folder, python implementation took about 7.5 seconds to run (Windows), breakdown was 100 ms reading and parsing file, 5000 ms getting kmers from file lines and histogramming, and then 2200 milliseconds sorting histogram contents.

Most of the kmers are unique in the file I tested.  There are also 'N' letters (!), not listed on wikipedia page, other sources imply these letters are when that DNA letter could not be read.

# c++ implementation
Plan is to use a hashmap / dictionary to map from kmer to count, but want to minimize sort time.  Can do this by using a priority queue to track only the topcount biggest values when outputting the most seen kmers.

Would expect the runtime of the various parts to be:

* File read: O(size) where size is the filesize
* Histogramming: O(uniques) where uniques is the total number of kmers (using a hashmap / dictionary)
* Topcount: O(uniques + topcount * (lg topcount) * (lg uniques)) where topcount is the number to output (probabilistic)
