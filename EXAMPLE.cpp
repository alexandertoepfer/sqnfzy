#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <vector>

#include <sqnmanip/sqn/fzy.hpp>

/**
 * readFasta
 *
 * @brief
 * Read complete file contents into a string.
 * @param fileName Name of the file in current working directory.
 * @return The file contents as single string.
 */
std::string
readFasta (const char* fileName)
{
    std::ifstream file (fileName);
    std::stringstream buffer;
    std::string::size_type index = 0;
    for (std::string line; getline (file, line);)
    {
        /* TODO: Fix Fasta EOF bug -> Ticket for Genbank format */
        index = line.find('\n', index);
        if (index != std::string::npos)
            line.erase(index);
        buffer << line;
    }
    return buffer.str ();
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
    sout << std::right << std::setw (5) << "[" + std::to_string (item._start) << " ... ";
    sout << item._sequence.toString ();
    sout << " ... " << std::left << std::setw (5) << std::to_string (item._end) + ")";
    return sout.str ();
}

/**
 * helper_center
 *
 * @brief
 * Center a string before passing it to print.
 * @param str The string to be centered.
 * @param width The total width used to calculate center.
 * @return The centered string.
 */
std::string
helper_center (std::string str, int width)
{
    std::string result;
    for (int i = 0; i < (width / 2) - (str.length() / 2); ++i)
        result += " ";
    return result += str;
}

int
main ()
{
    std::string sequence = readFasta ("data/AAV-CamKII-GCaMP6s-WPRE-SV40.fasta");
    sqn::Sequence<Dna5> genome (sequence), t7tag = "atggctagcatgactggtggacagcaaatgggt";
    
    sqn::FuzzyQuery<Dna5Sequence> analysis = { genome, t7tag };
    
    analysis.initializeScoreMatrix (sqn::continuityMatrix, 5);
    analysis.setItemParser (itemParse);
    
    for (auto match : analysis.search ())
    {
        std::cout << match.needle () << std::endl;
        std::cout << match.haystack () << std::endl;
        std::cout << helper_center (match.score (), match.haystack ().length ())
            << std::endl << std::endl;
    }
    
    return 0;
}
