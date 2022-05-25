#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include "SqnFzy.hpp"

/**
 * readFile
 *
 * @brief
 * Read complete file contents into a string.
 * @param fileName Name of the file in current working directory.
 * @return The file contents as single string.
 */
std::string
readFile (const char* fileName)
{
    std::ifstream t(fileName);
    std::stringstream buffer;
    buffer << t.rdbuf();
    return buffer.str();
}

/**
 * itemParse
 *
 * @brief
 * Parse a single item from a match as string.
 * @param item The item referring to a part of the sequence.
 * @return The formatted string.
 */
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

/**
 * center
 *
 * @brief
 * Center a string before passing it to print.
 * @param str The string to be centered.
 * @param width The total width used to calculate center.
 * @return The centered string.
 */
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
    
    return 0;
}
