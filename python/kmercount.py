import sys
from datetime import datetime

start = datetime.now()

if len(sys.argv) == 4:
    fastqFilename = sys.argv[1]
    k = int(sys.argv[2])
    topcount = int(sys.argv[3])
else:    
    fastqFilename = r"C:\code\fastq\fastq\example\ERR055763.filt.fastq"
    k = 30
    topcount = 25

verbose = True
linesProcessed = 0

with open(fastqFilename) as f:
    contents = f.read()

lines = [line.strip() for line in contents.split('\n') if line != '']

endfile = datetime.now()

countFromKmer = {}


for line in lines[1::4]:
    startIdx = 0
    while startIdx + k <= len(line):
        kmer = line[startIdx:(startIdx+k)]
        countFromKmer[kmer] = countFromKmer.get(kmer, 0) + 1
        startIdx += 1
    linesProcessed += 1

endhist = datetime.now()    
            
kmerCounts = list(countFromKmer.iteritems())
kmerCounts.sort(reverse = True, key = lambda kmerCount: kmerCount[1])

if len(kmerCounts) < topcount:
    raise Exception('fewer kmers than topcount')

end = datetime.now()

if verbose:
    uniqueKmers = len(countFromKmer.keys())
    totalKmers = sum(countFromKmer.values())
    print "=== Arguments ==="
    print "Filename: %s" % fastqFilename
    print "k: %d" % k
    print "top: %d" % topcount
    print 
    print "=== Detailed Information ==="
    print "DNA lines: %d" % linesProcessed
    print "Unique kmers: %d" % uniqueKmers
    print "Total kmers: %d" % totalKmers
    print "%% kmers unique: %3.1f" % (100.0 * uniqueKmers / totalKmers, )
    print "Time [ms]: %0.0f" % ((end - start).total_seconds() * 1000, )
    print "Time file [ms]: %0.0f" % ((endfile - start).total_seconds() * 1000, )
    print "Time hist [ms]: %0.0f" % ((endhist - endfile).total_seconds() * 1000, )
    print "Time sort [ms]: %0.0f" % ((end - endhist).total_seconds() * 1000, )
    print
    print "=== Top %d kmers ===" % topcount

for i in range(topcount):
    print '%s %d' % kmerCounts[i]
