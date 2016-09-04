// fastq.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <tuple>
#include <ctime>

void printTopKmers(size_t topcount, const std::unordered_map<std::string, int>& countFromKmer)
{
  std::priority_queue<std::tuple<int, std::string>, std::vector<std::tuple<int,  std::string>>, std::greater<std::tuple<int, std::string>>> pq;
  std::stack<std::tuple<int, std::string>> stack; // for reversing the p-queue
  
  std::tuple<int, std::string> next; // next item to be added to queue

  int count = 0;
  int bottomCount = 0;

  // first topcount iterations, expected work is (lg topcount)
  // topcount to numKeys(countFromKmer) iterations (on iteration n), expected work is 1 + topcount (lg topcount) / n
  // so topcount (lg topcount) (lg uniques) + uniques
  for (auto const& kvp : countFromKmer)
  {
    count = kvp.second;

    if ((pq.size() < topcount) || (bottomCount < count))
    {
      next = std::make_tuple(count, kvp.first);
      pq.push(next);
      if (pq.size() == 1 + topcount)
      {
        pq.pop();
        bottomCount = std::get<0>(pq.top());
      }
    }
  }

  //they're backwards, use a stack to put them forwards

  while(!pq.empty())
  {
    next = pq.top();
    stack.push(next);
    pq.pop();
  }

  while (!stack.empty())
  {
    next = stack.top();
    stack.pop();

    std::cout << std::get<1>(next) << " " << std::get<0>(next) << std::endl;
  }
}

void countKmers(size_t k, const std::vector<std::string>& lines, std::unordered_map<std::string, int>& countFromKmer)
{
  size_t lineSize = 0;
  std::string sub;
  for (auto const& line : lines)
  {
    lineSize = line.size();
    for (int i = 0 ; i + k <= lineSize ; ++i)
    {
      sub = line.substr(i, k);
      if (0 == countFromKmer.count(sub))
      {
        countFromKmer[sub] = 1;
      }
      else
      {
        countFromKmer[sub] += 1;
      }
    }
  }
}

int nonStrictParseFastq(const std::string& filename, std::vector<std::string>& lines)
{
  int linesProcessed = 0;
  std::ifstream inp(filename);
  int stanza = 0;
  for (std::string line; std::getline(inp, line) ; )
  {
    if (1 == stanza)
    {
      lines.push_back(line);
      linesProcessed += 1;
    }
    stanza = (stanza + 1) % 4;
  }
  inp.close();
  return linesProcessed;
}

int main(int argc, char **argv)
{
  if (4 != argc)
  {
    std::cout << "Usage: " << argv[0] << " filename k n" << std::endl;
    std::cout << "Example: " << argv[0] << " foo.fastq 25 30" << std::endl;
    exit(1);    
  }

  bool verbose = false;
  bool blockAtEnd = false;
  size_t k = atoi(argv[2]);
  size_t topcount = atoi(argv[3]);
  std::string fastqFilename(argv[1]);

  int linesProcessed = 0;
  clock_t begin;
  clock_t end;

  std::vector<std::string> lines;

  begin = clock();

  // read file in

  linesProcessed = nonStrictParseFastq( fastqFilename, lines );

  if (verbose)
  {
    std::cout << "read time " << double(clock() - begin) / (CLOCKS_PER_SEC / 1000) << std::endl;
    std::cout << "lines read: " << linesProcessed << std::endl;
    begin = clock();
  }

  size_t sizeGuess = lines.size() * (lines[0].length() - k);

  std::unordered_map<std::string, int> countFromKmer(sizeGuess);

  // count unique kmers

  countKmers(k, lines, countFromKmer);

  if (verbose)
  {
    std::cout << "hist time " << double(clock() - begin) / (CLOCKS_PER_SEC / 1000) << std::endl;
    begin = clock();
  }

  // grab top kmers

  printTopKmers(topcount, countFromKmer);

  if (verbose)
  {
    std::cout << "p-q time " << double(clock() - begin) / (CLOCKS_PER_SEC / 1000) << std::endl;
    begin = clock();
  }

  if (blockAtEnd)
  {
    getchar();
  }

  end = clock();

  return 0;
}

