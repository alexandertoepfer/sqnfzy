<p align="center">
  <img src="sqnfzy-logo-300.png" width="350">
</p>

# Analysis of DNA Sequences

Inspired by [University of Alberta, Paul Stothard's Sequence Manipulation Suite](https://www.bioinformatics.org/sms2/index.html) and [Seqan3](https://github.com/seqan/seqan3) a small bioinformatics library written in C++ to approximately find parts of DNA sequences that can be easily mutated into a useful restriction site. Currently this project is not open source.

# Usage

Initialization of sequences from string can be accomplished by assignment to a Sequence<Tp> instance with the given DNA or RNA sequence type.
```c
sqn::Sequence<Dna5> randGenome = "CTTTACAGGCCCCGGTTTCT";
sqn::Sequence<Dna5> enzymeEagI = "CGGCCG";
```

Initialization of the fuzzy query with genome (haystack) and enzyme (needle) according to their DNA or RNA sequence type.
```c
sqn::FuzzyQuery<Dna5Sequence> query = { randGenome, enzymeEagI };
```

Configuration of the query matrix used to calculate scores for sequence parts. 
```c
sqn::ScoreMatrix scores = { /* Match */ 1, /* Mismatch */ -1, /* Gap */ 2 };
query.initializeScoreMatrix(scores, /* Amount of matches */ 3);
```

Execution of the query and collection of results.
```c
for (sqn::Match<Dna5Sequence>& m : query.search ())
{
  std::string needle = m.needle ();
  std::string haystack = m.haystack ();
  //std::string score = m.score ();
}
```

> **Note:** For custom string formats define std::string itemParse (sqn::Item<Tp>& item) 

# Example
  
Code
