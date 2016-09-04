// fastq.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h"
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

void printTopKmers(unsigned int topcount, const std::unordered_map<std::string, int>& countFromKmer)
{
	auto cmp = [](std::tuple<std::string, int> left, std::tuple<std::string, int> right) { return std::get<1>(left) > std::get<1>(right) ; };
	std::priority_queue<std::tuple<std::string, int>, std::vector<std::tuple<std::string, int>>, decltype(cmp)> pq;
	std::stack<std::tuple<std::string, int>> stack; // for reversing the p-queue
	
	std::tuple<std::string, int> next; // next item to be added to queue

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
			next = std::make_tuple(kvp.first, count);
			pq.push(next);
			if (pq.size() == 1 + topcount)
			{
				pq.pop();
				bottomCount = std::get<1>(pq.top());
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

		std::cout << std::get<0>(next) << " " << std::get<1>(next) << std::endl;
	}
}

void countKmers(unsigned int k, const std::vector<std::string>& lines, std::unordered_map<std::string, int>& countFromKmer)
{
	unsigned int lineSize = 0;
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

int _tmain(int argc, _TCHAR* argv[])
{
	bool verbose = true;
	bool blockAtEnd = false;
	unsigned int k = 30;
	unsigned int topcount = 25;
	std::string fastqFilename("C:\\Users\\mfishburn\\Downloads\\ERR055763.filt.fastq");

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

	unsigned int sizeGuess = lines.size() * (lines[0].length() - k);

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

