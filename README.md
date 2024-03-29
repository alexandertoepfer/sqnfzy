# sqnfzy - Fuzzy Analysis of DNA Sequences

[![Codacy Badge](https://app.codacy.com/project/badge/Grade/fee733701e2943a19d93515d0104c2ef)](https://www.codacy.com/gh/alexandertoepfer/sqnfzy/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=alexandertoepfer/sqnfzy&amp;utm_campaign=Badge_Grade) ![GitHub](https://img.shields.io/github/license/alexandertoepfer/sqnfzy?style=flat)

Inspired by [University of Alberta, Paul Stothard's Work](https://www.bioinformatics.org/sms2/index.html) and [Seqan3](https://github.com/seqan/seqan3) a small bioinformatics library written in C++ to approximately find parts of DNA sequences that can be easily mutated into a useful restriction site.

> **Note:** In case you want to make use of this algorithm, I highly recommend using a faster implementation of levenstein at [RapidFuzz](https://github.com/maxbachmann/RapidFuzz)

# Usage

Initialization of sequences from string can be accomplished by assignment to a sequence instance with the given DNA or RNA sequence type.
```c
sqn::Sequence<Dna5> randGenome = "CTTTACAGGCCCCGGTTTCT";
sqn::Sequence<Dna5> enzymeEagI = "CGGCCG";
```

Initialization of the fuzzy query with genome (haystack) and enzyme (needle) according to their DNA or RNA sequence type.
```c
sqn::FuzzyQuery<Dna5Sequence> query = {randGenome, enzymeEagI};
```

Configuration of the values used to calculate scores for sequence parts can be done by setting a score matrix that dictates rewards and penalties.
```c
sqn::ScoreMatrix scores = {/*match*/1, /*mismatch*/-1, /*gap*/2 };
query.initializeScoreMatrix(scores, /*amount of matches*/3);
```

For execution of the query and collection of results it is generally recommended to wrap everything into a loop, as std::list<Match<Tp>>& **search** () will return a list of all approximate matches.
```c
for(sqn::Match<Dna5Sequence>& match : query.search()) {
  std::string needle = match.needle();
  std::string haystack = match.haystack();
}
```

The score of the match can be obtained as well, which is shown in the complete example but shall be left out here. 
> **Note:** For custom string formats of matches define std::string **itemParse** (sqn::Item<Tp>& item) and pass it as function pointer to the query.

# Example

Using the functionality provided by the library to find possible common plasmid features like T7 in adeno-associated virus sequences allowing for errors and gaps caused by mutations. Take a look at <code>example.cpp</code> for the complete code. Sequence used in this example [Addgene #107790-AAV9](https://www.addgene.org/browse/sequence/204876/)
```c
#include <sqnmanip/sqn/fzy.hpp>
  
std::string itemParse (sqn::Item<Dna5Sequence>& item) {
  ...
}

int main() {
  sqn::Sequence<Dna5> genome(readFasta("AAV-CamKII-GCaMP6s-WPRE-SV40.fasta"));
  sqn::Sequence<Dna5> t7tag = "atggctagcatgactggtggacagcaaatgggt";
  sqn::FuzzyQuery<Dna5Sequence> analysis = {genome, t7tag};
  
  analysis.initializeScoreMatrix(sqn::continuityMatrix, 3);
  analysis.setItemParser(itemParse);
  
  for(sqn::Match<Dna5Sequence>& match : analysis.search()) {
    ...
  }
  return 0;
}
```

Using the functionality provided by the library to auto correct text based on a dictionary of known words, currently not officially supported but another interesting use case for this algorithm.
```c
#include <sqnmanip/sqn/fzy.hpp>

int main() {
  sqn::Sequence<char> words = "student,summer,school,system,sample";
  sqn::Sequence<char> word = "stwdnt";
  
  sqn::FuzzyQuery<sqn::Sequence<char>> autocorrect = {words, word};
  autocorrect.initializeScoreMatrix(sqn::standardMatrix, 1);
  
  for (sqn::Match<sqn::Sequence<char>>& match : autocorrect.search()) {
    ...
  }
  return 0;
}
```
  
# Results
From two different score matrices besides the obvious first match multiple other approximate matches could be found even though they are scoring below half of the perfect score.
```c
                       (DISPARITY MATRIX)
              Gaps have less penalty than mismatches
    
         [1 ... ATGGCTAGCATGACTGGTGGACAGCAAATGGGT ... 33)
      [1509 ... ATGGCTAGCATGACTGGTGGACAGCAAATGGGT ... 1541)
                            Score 33 (Perfect)
    
         [4 ... GCT-AGCATGACTGGTG-GACA-GCA-AATGG ... 31)
       [689 ... GCTAAG-GTGGC-GGTGTGATATGCACAATGG ... 718)
                            Score 14
    
       [1 ... ATGGCTAGCATGACTGGTGGA--CA-GCAA--ATGG ... 31)
    [5579 ... AT--CTA-CACGAC-GG-GGAGTCAGGCAACTATGG ... 5609)
                            Score 14
     

                       (CONTINUITY MATRIX)
           Mismatches have low penalty, gaps high penalty
    
         [1 ... ATGGCTAGCATGACTGGTGGACAGCAAATGGGT ... 33)
      [1509 ... ATGGCTAGCATGACTGGTGGACAGCAAATGGGT ... 1541)
                            Score 33 (Perfect)
    
         [1 ... ATGGCTAGCATGACTGGTGGACAGCAAATGGG ... 32)
      [1323 ... AAGGCTCGCGAGGCT-GTGAGCAGCCACAGTG ... 1353)
                            Score 18
    
         [1 ... ATGGCTAGCATGACTGGTGGACAGCAAATGGG ... 32)
      [3170 ... CTGCCTTGCCCG-CTGCTGGACAGGGGCTCGG ... 3200)
                            Score 18
     
```

Results for auto correction of text with a given dictionary as follows.
```c
                       (STANDARD MATRIX)
           Mismatches and gaps have similar penalties
                   
                      stwd-nt from 1 to 6
                      student from 1 to 7
                            Score 3
```
