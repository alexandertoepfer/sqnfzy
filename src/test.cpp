#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include "SqnFzy.hpp"

std::string
readFile (const char* fileName)
{
    std::ifstream t(fileName);
    std::stringstream buffer;
    buffer << t.rdbuf();
    return buffer.str();
}

template<typename Tp>
std::string
itemParse (sqn::Item<Tp>& item)
{
    std::stringstream sout;
    sout << std::right << std::setw (5) << "[" + std::to_string(item._start) << " ... ";
    sout << item._sequence.toString ();
    sout << " ... " << std::left << std::setw (5) << std::to_string(item._end) + ")";
    return sout.str ();
}

std::string
center (std::string str, int width)
{
    std::string result;
    for (int i = 0; i < (width / 2) - (str.length() / 2); ++i)
        result += " ";
    return result += str;
}

int
main ()
{
    sqn::Sequence<Dna5> genome(readFile("AAV-CamKII-GCaMP6s-WPRE-SV40.fasta"));
    sqn::Sequence<Dna5> t7tag = "atggctagcatgactggtggacagcaaatgggt";
    sqn::FuzzyQuery<Dna5Sequence> analysis = { genome, t7tag };
    analysis.initializeScoreMatrix(sqn::continuityMatrix, 3);
    analysis.setItemParser(itemParse);
    for (sqn::Match<Dna5Sequence>& m : analysis.search ())
    {
        std::cout << m.needle () << std::endl;
        std::cout << m.haystack () << std::endl;
        std::cout << center (m.score (), m.haystack ().length ())
            << std::endl << std::endl;
    }

    /*                    (DISPARITY MATRIX)
     *           Gaps have less penalty than mismatches
     *
     *     [1 ... ATGGCTAGCATGACTGGTGGACAGCAAATGGGT ... 33)
     *  [1509 ... ATGGCTAGCATGACTGGTGGACAGCAAATGGGT ... 1541)
     *                        Score 33 (Perfect)
     *
     *     [4 ... GCT-AGCATGACTGGTG-GACA-GCA-AATGG ... 31)
     *   [689 ... GCTAAG-GTGGC-GGTGTGATATGCACAATGG ... 718)
     *                        Score 14
     *
     *     [1 ... ATGGCTAGCATGACTGGTGGA--CA-GCAA--ATGG ... 31)
     *  [5579 ... AT--CTA-CACGAC-GG-GGAGTCAGGCAACTATGG ... 5609)
     *                          Score 14
     *
     *     [2 ... TGGCTAGCATGACTGGTGGACAGCAAATGGGT ... 33)
     *   [603 ... T-GCT-GC-TCAGTGGT-GACAG-ATAGGGGT ... 629)
     *                        Score 14
     *
     *     [1 ... AT-G-GCT--A--GCAT-GAC-T---GGTGGACAGCAA-AT ... 29)
     *  [2316 ... ATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACAT ... 2356)
     *                            Score 13
     */

    /*                   (CONTINUITY MATRIX)
     *       Mismatches have low penalty, gaps high penalty
     *
     *     [1 ... ATGGCTAGCATGACTGGTGGACAGCAAATGGGT ... 33)
     *  [1509 ... ATGGCTAGCATGACTGGTGGACAGCAAATGGGT ... 1541)
     *                        Score 33 (Perfect)
     *
     *     [1 ... ATGGCTAGCATGACTGGTGGACAGCAAATGGG ... 32)
     *  [1323 ... AAGGCTCGCGAGGCT-GTGAGCAGCCACAGTG ... 1353)
     *                        Score 18
     *
     *     [1 ... ATGGCTAGCATGACTGGTGGACAGCAAATGGG ... 32)
     *  [3170 ... CTGCCTTGCCCG-CTGCTGGACAGGGGCTCGG ... 3200)
     *                        Score 18
     */

    return 0;
}
