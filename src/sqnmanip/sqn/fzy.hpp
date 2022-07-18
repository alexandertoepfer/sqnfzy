#include <algorithm>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <fstream>

#ifndef SQN_FZY_HPP
#define SQN_FZY_HPP

/*
 * MIT License
 * 
 * Copyright (c) 2022 Alexander Töpfer alexander_toepfer@gmx.de
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*
 * GNU General Public License
 * 
 * Copyright (C) 2020 Paul Stothard stothard@ualberta.ca
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
 
namespace sqn
{
    template <class T, class U>
    concept DerivedFrom = std::is_base_of<U, T>::value;

    template <class U, class T>
    concept AssignableTo = std::is_assignable<T&, U>::value;

    /**
     * SequencePointerBase
     *
     * @brief
     * Base class for navigation between sequence nodes of any char type elements.
     */
    struct SequencePointerBase
    {
        SequencePointerBase* _next;
        SequencePointerBase* _prev;
    };

    /**
     * SequenceNodeBase
     *
     * @brief
     * Base class to contruct sequence node from any char type element.
     */
    template <typename Tp>
    struct SequenceNodeBase : public SequencePointerBase
    {
        Tp _data;
    };

    /**
     * DnaSpecification
     *
     * @brief
     * Base class for DNA alphabet types.
     */
    class DnaSpecification
    {
      public:
        virtual std::map<size_t, unsigned char>& nucleotides () = 0;
    };

    /**
     * SequenceIterator
     *
     * @brief
     * Base class template has the ability to iterate through the elements of sequences.
     */
    template <typename Tp>
    class SequenceIterator
    {
      public:
        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type = std::ptrdiff_t;

        SequenceIterator (SequenceNodeBase<Tp>* position)
            : _current (position) {};

        /**
         * operator==, operator!=
         *
         * @brief
         * Comparing whether iterators are pointing to the same element.
         * @param other The second iterator for comparison.
         */
        bool
        operator== (SequenceIterator<Tp> other) const
        {
            return this->_current == other._current;
        }
        bool
        operator!= (SequenceIterator<Tp> other) const
        {
            return this->_current != other._current;
        }

        /**
         * operator+, operator-
         *
         * @brief
         * Iterate an iterator forward by addition and backward by subtraction.
         * @param summand, subtrahend The amount of forward backward operations.
         */
        SequenceIterator<Tp>
        operator+ (unsigned int summand)
        {
            SequenceIterator<Tp> tmp (this->_current);
            while (summand > 0)
            {
                tmp._current
                    = static_cast<SequenceNodeBase<Tp>*> (tmp._current->_next);
                --summand;
            }
            return tmp;
        }
        SequenceIterator<Tp>
        operator- (unsigned int subtrahend)
        {
            SequenceIterator<Tp> tmp (this->_current);
            while (subtrahend > 0)
            {
                tmp._current
                    = static_cast<SequenceNodeBase<Tp>*> (tmp._current->_prev);
                --subtrahend;
            }
            return tmp;
        }

        /**
         * operator++, operator--
         *
         * @brief
         * Iterate an iterator forward or backward, this would be equal to either
         * operator+(1) or operator-(1) depending on prefix or postfix.
         */
        SequenceIterator<Tp>&
        operator++ ()
        {
            this->_current
                = static_cast<SequenceNodeBase<Tp>*> (this->_current->_next);
            return *this;
        }
        SequenceIterator<Tp>&
        operator-- ()
        {
            this->_current
                = static_cast<SequenceNodeBase<Tp>*> (this->_current->_prev);
            return *this;
        }
        SequenceIterator<Tp>
        operator++ (int)
        {
            SequenceIterator<Tp> tmp (this->_current);
            this->operator++ ();
            return tmp;
        }
        SequenceIterator<Tp>
        operator-- (int)
        {
            SequenceIterator<Tp> tmp (this->_current);
            this->operator-- ();
            return tmp;
        }

        /**
         * operator*, operator->
         *
         * @brief
         * Dereference the current iterator and access the node element.
         */
        Tp&
        operator* () const
        {
            return this->_current->_data;
        }
        Tp*
        operator-> ()
        {
            return this->_current;
        }

      private:
        SequenceNodeBase<Tp>* _current;
    };

    /**
     * Sequence
     *
     * @brief
     * Class template that represents sequences of DNA alphabet types.
     */
    template <AssignableTo<char> Tp>
    class Sequence
    {
      public:
        Sequence ()
            : _begin (nullptr), _end (nullptr), _size (0) {};
        Sequence (const char* other)
        {
            this->_size = 0;
            this->_begin = nullptr;
            this->_end = nullptr;
            while (other[this->_size] != '\0')
            {
                SequenceNodeBase<Tp>* node = new SequenceNodeBase<Tp> ();
                node->_data = other[this->_size];
                this->append (node);
                ++this->_size;
            }
        }
        Sequence (std::string other)
        {
            this->_size = 0;
            this->_begin = nullptr;
            this->_end = nullptr;
            for (std::string::size_type i = 0; i < other.size(); i++)
            {
                SequenceNodeBase<Tp>* node = new SequenceNodeBase<Tp> ();
                node->_data = other[i];
                this->append (node);
                ++this->_size;
            }
        }
        Sequence (const Sequence<Tp>& other)
        {
            this->_begin = other._begin;
            this->_end = other._end;
            this->_size = other._size;
        }

        /**
         * operator+=
         *
         * @brief
         * Append to the current sequence.
         * @param summand The element or sequence to append.
         */
        Sequence<Tp>&
        operator+= (const char* summand)
        {
            size_t index = 0;
            while (summand[index] != '\0')
            {
                SequenceNodeBase<Tp>* node = new SequenceNodeBase<Tp> ();
                node->_data = summand[index];
                this->append (node);
                ++this->_size;
                ++index;
            }
            return *this;
        }
        template <typename Ts>
        Sequence<Tp>&
        operator+= (Ts summand)
        {
            SequenceNodeBase<Tp>* node = new SequenceNodeBase<Tp> ();
            node->_data = summand;
            this->append (node);
            ++this->_size;
            return *this;
        }

        /**
         * operator[]
         *
         * @brief
         * Access the current sequence element at a specified index.
         * @param index The position in the sequence.
         */
        Tp&
        operator[] (size_t index)
        {
            SequenceNodeBase<Tp>* node = this->_begin;
            for (size_t i = 0; i < index; ++i)
                node = static_cast<SequenceNodeBase<Tp>*> (node->_next);
            return node->_data;
        }

        /**
         * begin(), end()
         *
         * @brief
         * Retrieve an iterator pointing to the begin or end of the sequence.
         * @return The SequenceIterator with that position.
         */
        SequenceIterator<Tp>
        begin ()
        {
            return SequenceIterator<Tp> (this->_begin);
        }
        SequenceIterator<Tp>
        end ()
        {
            return SequenceIterator<Tp> (this->_end);
        }
        size_t
        length ()
        {
            return this->_size;
        }
        void
        reverse ()
        {
            auto first = this->begin ();
            auto last = this->end ();

            if (first == last - 1 || first + 1 == last)
                return std::swap (*first, *last);

            while (first != last && first != last - 1)
                std::swap (*first++, *last--);

            std::swap (*first, *last);
        }
        std::string
        toString ()
        {
            std::string result;
            SequenceNodeBase<Tp>* node = this->_begin;
            while (node != nullptr)
            {
                result += (char)node->_data;
                node = static_cast<SequenceNodeBase<Tp>*> (node->_next);
            }
            return result;
        }

      private:
        void
        append (SequenceNodeBase<Tp> *node)
        {
            if (this->_end != nullptr)
            {
                this->_end->_next = node;
                node->_prev = this->_end;
                this->_end = node;
            }
            else
            {
                this->_begin = node;
                this->_end = node;
            }
        }
        size_t _size;
        SequenceNodeBase<Tp>* _begin;
        SequenceNodeBase<Tp>* _end;
    };

    /**
     * Dna5Specification
     *
     * @brief
     * The five letter DNA alphabet of A,C,G,T, the unknown character N and a gap character.
     */
    class Dna5Specification : DnaSpecification
    {
      public:
        std::map<size_t, unsigned char>&
        nucleotides ()
        {
            return this->_dna5nucleotides;
        }

      private:
        std::map<size_t, unsigned char> _dna5nucleotides
            = { { 65, 'A' }, { 67, 'C' }, {  71, 'G' }, {  84, 'T' }, { 45, '-' },
                { 97, 'A' }, { 99, 'C' }, { 103, 'G' }, { 116, 'T' } };
    };

    /**
     * SimpleType
     *
     * @brief
     * Class template that stores the values of an instance according to its specification.
     */
    template <typename TValue, DerivedFrom<DnaSpecification> TSpec>
    class SimpleType
    {
      public:
        TValue
        getNucleotides (const size_t& key,
                        const TValue& value)
        {
            typename std::map<size_t, TValue>::const_iterator it
                = this->_specification.nucleotides ().find (key);
            if (it == this->_specification.nucleotides ().end ())
                return value;
            else
                return it->second;
        }

        /**
         * operator=
         *
         * @brief
         * Assign a value with specification restriction.
         */
        template <typename T>
        inline SimpleType&
        operator= (T const& other)
        {
            this->_value = this->getNucleotides ((size_t)other, static_cast<TValue> ('N'));
            return *this;
        }

        operator char () const { return (char)this->_value; }

      private:
        TValue _value;
        TSpec _specification;
    };

    /**
     * Item
     *
     * @brief
     * Class template that stores part of a sequence with start and end.
     */
    template <typename Tp>
    class Item
    {
      public:
        Item () = default;
        Item (Tp sequence, size_t start, size_t end)
            : _sequence (sequence), _start (start), _end (end)
        {
        }
        Item (const Item &other)
            : _sequence (other._sequence), _start (other._start),
              _end (other._end)
        {
        }

        Tp _sequence;
        size_t _start;
        size_t _end;
    };

    /**
     * Match
     *
     * @brief
     * Class template that stores a sequence match with haystack and needle.
     */
    template <typename Tp>
    class Match
    {
      public:
        Match () = default;
        Match (Tp sequenceNeedle, Tp sequenceHaystack, int score,
               size_t startNeedle, size_t endNeedle,
               size_t startHaystack, size_t endHaystack,
               std::string (*parser)(Item<Tp>&))
        {
            this->_haystack
                = Item<Tp> (sequenceHaystack, startHaystack, endHaystack);
            this->_needle
                = Item<Tp> (sequenceNeedle, startNeedle, endNeedle);
            this->_score = score;
            this->_parser = parser;
        }

        /**
         * haystack()
         *
         * @brief
         * Get part of the haystack sequence that was matched with the needle.
         * @return String representation of the match.
         */
        std::string
        haystack ()
        {
            if (this->_parser)
            {
                return this->_parser(this->_haystack);
            }
            else
            {
                return this->_haystack._sequence.toString () + ' ' +
                    "from " + std::to_string(this->_haystack._start) + " to "
                    + std::to_string(this->_haystack._end);
            }
        }

        /**
         * needle()
         *
         * @brief
         * Get part of the needle sequence that was matched with the haystack.
         * @return String representation of the match.
         */
        std::string
        needle ()
        {
            if (this->_parser)
            {
                return this->_parser(this->_needle);
            }
            else
            {
                return this->_needle._sequence.toString () + ' ' +
                    "from " + std::to_string(this->_needle._start) + " to "
                    + std::to_string(this->_needle._end);
            }
        }

        /**
         * score()
         *
         * @brief
         * Get the calculated score for the match.
         */
        std::string
        score()
        {
            std::string res = "Score ";
            res += std::to_string(this->_score);
            return res;
        }

        int _score;
        Item<Tp> _haystack;
        Item<Tp> _needle;
        std::string (*_parser)(Item<Tp>&);
    };

    /**
     * ScoreMatrix
     *
     * @brief
     * Class that describes how fuzzy scores are calculated when comparing two sequence characters.
     */
    class ScoreMatrix
    {
      public:
        ScoreMatrix () = default;
        ScoreMatrix (const ScoreMatrix& other)
        {
            this->_match = other._match;
            this->_mismatch = other._mismatch;
            this->_gap = other._gap;
        }
        ScoreMatrix (int match, int mismatch, int gap)
        {
            this->_match = match;
            this->_mismatch = mismatch;
            this->_gap = gap;
        }

        /**
         * getScore()
         *
         * @brief
         * Calculate the score for two characters from the haystack and needle.
         * @param r1 Character from haystack.
         * @param r2 Character from needle.
         */
        int
        getScore (char r1, char r2)
        {
            if (!(tolower (r1) - tolower (r2)))
            {
                return this->_match;
            }
            else
            {
                return this->_mismatch;
            }
        }

        int _match;
        int _mismatch;
        int _gap;
    };

    /**
     * Node
     *
     * @brief
     * Class that describes one node unit to span the haystack*needle matrix,
     * it stores the calculated fuzzy score and position.
     */
    class Node
    {
      public:
        Node () : _value (Node::_undefined),
                  _tracebackI (Node::_undefined),
                  _tracebackJ (Node::_undefined) {}

        int _value;
        int _tracebackI;
        int _tracebackJ;
        bool _alreadyMatched = false;
        static const int _undefined = std::numeric_limits<int>::max();
    };

    /**
     * FuzzyQuery
     *
     * @brief
     * Class template fuzzy search algorithm calculates a fuzzy score for all nodes in matrix
     * and returns matching parts of the haystack sequence according to their score.
     */
    template <typename Tp>
    class FuzzyQuery
    {
      public:
        FuzzyQuery(Tp haystackSequence, Tp needleSequence)
        {
            this->initializeMatrix(haystackSequence, needleSequence);
            this->_parser = nullptr;
        }

        /**
         * setItemParser()
         *
         * @brief
         * Set parser for a single item from a match.
         * @param parser The function used for parsing part of a sequence with start and end.
         */
        void
        setItemParser(std::string (*parser)(Item<Tp>&))
        {
            this->_parser = parser;
        }

        /**
         * initializeScoreMatrix()
         *
         * @brief
         * Set rewards and penalties for score calculation as well as amounts.
         * @param scoreSet The values used to calculate scores for sequence parts.
         * @param amount The amount of matches to be retrieved.
         */
        void
        initializeScoreMatrix (ScoreMatrix scoreSet, int amount)
        {
            this->_scoreSet = scoreSet;
            this->_amount = amount;
        }

        /**
         * initializeMatrix()
         *
         * @brief
         * Span the matrix of haystack*needle dimensions and set undefined states.
         * @param haystackSequence The complete genome to be searched.
         * @param needleSequence The needle sequence to be approximately searched.
         */
        void
        initializeMatrix (Tp haystackSequence, Tp needleSequence)
        {
            this->_haystackSequence = haystackSequence;
            this->_needleSequence = needleSequence;
            this->_score = 0;

            std::vector<std::vector<Node>> nodes (
                needleSequence.length () + 1,
                std::vector<Node> (haystackSequence.length () + 1, Node ()));

            nodes[0][0]._value = 0;

            for (int i = 1; i < nodes.size (); i++)
            {
                nodes[i][0]._value = nodes[i - 1][0]._value;
                nodes[i][0]._tracebackI = i - 1;
                nodes[i][0]._tracebackJ = 0;
            }

            for (int j = 1; j < nodes[0].size (); j++)
            {
                nodes[0][j]._value = nodes[0][j - 1]._value;
                nodes[0][j]._tracebackI = 0;
                nodes[0][j]._tracebackJ = j - 1;
            }

            this->_nodes = nodes;
        }

        /**
         * updateMatrix()
         *
         * @brief
         * Update all nodes in the matrix.
         */
        void
        updateMatrix ()
        {
            for (int i = 1; i < this->_nodes.size (); i++)
            {
                for (int j = 1; j < this->_nodes[0].size (); j++)
                {
                    int a, b, c;

                    if (this->_nodes[i][j]._alreadyMatched)
                    {
                        a = 0;
                        b = 0;
                        c = 0;
                    }
                    else if (i == (this->_nodes.size () - 1)
                        && j == (this->_nodes[0].size () - 1))
                    {
                        a = this->_nodes[i - 1][j]._value;
                        b = this->_nodes[i][j - 1]._value;
                        c = this->_nodes[i - 1][j - 1]._value
                            + this->_scoreSet.getScore (this->_needleSequence[i - 1],
                            this->_haystackSequence[j - 1]);
                    }
                    else if (i == (this->_nodes.size () - 1))
                    {
                        a = this->_nodes[i - 1][j]._value
                            - this->_scoreSet._gap;
                        b = this->_nodes[i][j - 1]._value;
                        c = this->_nodes[i - 1][j - 1]._value
                            + this->_scoreSet.getScore (this->_needleSequence[i - 1],
                            this->_haystackSequence[j - 1]);
                    }
                    else if (j == (this->_nodes[0].size () - 1))
                    {
                        a = this->_nodes[i - 1][j]._value;
                        b = this->_nodes[i][j - 1]._value
                            - this->_scoreSet._gap;
                        c = this->_nodes[i - 1][j - 1]._value
                            + this->_scoreSet.getScore (this->_needleSequence[i - 1],
                            this->_haystackSequence[j - 1]);
                    }
                    else
                    {
                        a = this->_nodes[i - 1][j]._value
                            - this->_scoreSet._gap;
                        b = this->_nodes[i][j - 1]._value
                            - this->_scoreSet._gap;
                        c = this->_nodes[i - 1][j - 1]._value
                            + this->_scoreSet.getScore (this->_needleSequence[i - 1],
                            this->_haystackSequence[j - 1]);
                    }

                    if ((a > b) && (a > c))
                    {
                        this->_nodes[i][j]._value = a;
                        this->_nodes[i][j]._tracebackI = i - 1;
                        this->_nodes[i][j]._tracebackJ = j;
                    }
                    else if ((b > c) && (b > a))
                    {
                        this->_nodes[i][j]._value = b;
                        this->_nodes[i][j]._tracebackI = i;
                        this->_nodes[i][j]._tracebackJ = j - 1;
                    }
                    else
                    {
                        this->_nodes[i][j]._value = c;
                        this->_nodes[i][j]._tracebackI = i - 1;
                        this->_nodes[i][j]._tracebackJ = j - 1;
                    }
                    if (this->_nodes[i][j]._value < 0)
                    {
                        this->_nodes[i][j]._value = 0;
                        this->_nodes[i][j]._tracebackI = Node::_undefined;
                        this->_nodes[i][j]._tracebackJ = Node::_undefined;
                    }

                }
            }

            this->_score
                = this->_nodes[this->_nodes.size () - 1][this->_nodes[0].size () - 1]
                ._value;
        }

        /**
         * search()
         *
         * @brief
         * Get aligned haystack sequence and needle sequence matches according to amount.
         */
        std::list<Match<Tp>>&
        search ()
        {
            int hitCount = 0;
            while (hitCount < this->_amount)
            {
                this->updateMatrix ();
                int maxNodeValue = 0;
                int maxNodeI = Node::_undefined;
                int maxNodeJ = Node::_undefined;

                for (int i = 1; i < this->_nodes.size (); i++)
                {
                    for (int j = 1; j < this->_nodes[0].size (); j++)
                    {
                        if (this->_nodes[i][j]._value > maxNodeValue)
                        {
                            maxNodeValue = this->_nodes[i][j]._value;
                            maxNodeI = i;
                            maxNodeJ = j;
                        }
                    }
                }

                if (maxNodeValue == 0)
                {
                    break;
                }

                int currentI = maxNodeI;
                int currentJ = maxNodeJ;
                Node &currentNode = this->_nodes[maxNodeI][maxNodeJ];

                Tp alignedNeedle;
                Tp alignedHaystack;
                int score = currentNode._value;
                int endNeedle = Node::_undefined;
                int endHaystack = Node::_undefined;

                while (currentNode._tracebackI != Node::_undefined
                    && currentNode._tracebackJ != Node::_undefined)
                {
                    if (currentI == 0 || currentJ == 0)
                    {
                        break;
                    }

                    if (currentNode._tracebackI == currentI - 1
                        && currentNode._tracebackJ == currentJ - 1)
                    {
                        if (endNeedle == Node::_undefined)
                        {
                            endNeedle = currentI;
                            endHaystack = currentJ;
                        }
                        alignedNeedle += this->_needleSequence[currentI - 1];
                        alignedHaystack += this->_haystackSequence[currentJ - 1];
                    }
                    else if (currentNode._tracebackJ == currentJ - 1)
                    {
                        if (endNeedle != Node::_undefined)
                        {
                            alignedNeedle += "-";
                            alignedHaystack += this->_haystackSequence[currentJ - 1];
                        }
                    }
                    else
                    {
                        if (endNeedle != Node::_undefined)
                        {
                            alignedNeedle += this->_needleSequence[currentI - 1];
                            alignedHaystack += "-";
                        }
                    }

                    this->_nodes[currentI][currentJ]._value = 0;
                    this->_nodes[currentI][currentJ]._alreadyMatched = true;

                    currentI = currentNode._tracebackI;
                    currentJ = currentNode._tracebackJ;

                    currentNode = this->_nodes[currentNode._tracebackI][currentNode._tracebackJ];
                }

                alignedNeedle.reverse ();
                alignedHaystack.reverse ();

                Match<Tp> match (alignedNeedle, alignedHaystack, score,
                    currentI + 1, endNeedle, currentJ + 1, endHaystack,
                    this->_parser != nullptr ? this->_parser : nullptr);
                this->_matches.push_back (match);
                hitCount++;
            }
            return this->_matches;
        }

        private:
          Tp _needleSequence;
          Tp _haystackSequence;
          int _score;
          int _amount;
          ScoreMatrix _scoreSet;
          std::vector<std::vector<Node>> _nodes;
          std::list<Match<Tp>> _matches;
          std::string (*_parser)(Item<Tp>&);
    };

    ScoreMatrix continuityMatrix { 1, 0, 2 }, disparityMatrix { 1, -1, 0 }, standardMatrix { 1, -1, 1 };
}

typedef sqn::SimpleType<unsigned char, sqn::Dna5Specification> Dna5;
typedef sqn::Sequence<Dna5> Dna5Sequence;

#endif