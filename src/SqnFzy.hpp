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
        operator== (SequenceIterator<Tp> other) const;

        bool
        operator!= (SequenceIterator<Tp> other) const;

        /**
         * operator+, operator-
         *
         * @brief
         * Iterate an iterator forward by addition and backward by subtraction.
         * @param summand, subtrahend The amount of forward backward operations.
         */
        SequenceIterator<Tp>
        operator+ (unsigned int summand);

        SequenceIterator<Tp>
        operator- (unsigned int subtrahend);

        /**
         * operator++, operator--
         *
         * @brief
         * Iterate an iterator forward or backward, this would be equal to either
         * operator+(1) or operator-(1) depending on prefix or postfix.
         */
        SequenceIterator<Tp>&
        operator++ ();

        SequenceIterator<Tp>&
        operator-- ();

        SequenceIterator<Tp>
        operator++ (int);

        SequenceIterator<Tp>
        operator-- (int);

        /**
         * operator*, operator->
         *
         * @brief
         * Dereference the current iterator and access the node element.
         */
        Tp&
        operator* () const;

        Tp*
        operator-> ();

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
            : _begin (nullptr), _end (nullptr), _size (0) {}
        
        Sequence (const Sequence<Tp>& other);
        
        Sequence (const char* other);

        Sequence (std::string other);

        /**
         * operator+=
         *
         * @brief
         * Append to the current sequence.
         * @param summand The element or sequence to append.
         */
        Sequence<Tp>&
        operator+= (const char* summand);

        template <typename Ts>
        Sequence<Tp>&
        operator+= (Ts summand);

        /**
         * operator[]
         *
         * @brief
         * Access the current sequence element at a specified index.
         * @param index The position in the sequence.
         */
        Tp&
        operator[] (size_t index);

        /**
         * begin(), end()
         *
         * @brief
         * Retrieve an iterator pointing to the begin or end of the sequence.
         * @return The SequenceIterator with that position.
         */
        SequenceIterator<Tp>
        begin ();

        SequenceIterator<Tp>
        end ();

        size_t
        length ();

        void
        reverse ();

        std::string
        toString ();

      private:
        void
        append (SequenceNodeBase<Tp> *node);

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
        nucleotides ();

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
                        const TValue& value);

        /**
         * operator=
         *
         * @brief
         * Assign a value with specification restriction.
         */
        template <typename T>
        inline SimpleType&
        operator= (T const& other);


        operator char () const;

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
            : _sequence (sequence), _start (start), _end (end) {}
        Item (const Item &other)
            : _sequence (other._sequence), _start (other._start),
              _end (other._end) {}

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
               std::string (*parser)(Item<Tp>&));

        /**
         * haystack()
         *
         * @brief
         * Get part of the haystack sequence that was matched with the needle.
         * @return String representation of the match.
         */
        std::string
        haystack ();

        /**
         * needle()
         *
         * @brief
         * Get part of the needle sequence that was matched with the haystack.
         * @return String representation of the match.
         */
        std::string
        needle ();

        /**
         * score()
         *
         * @brief
         * Get the calculated score for the match.
         */
        std::string
        score();

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
        ScoreMatrix (const ScoreMatrix& other);

        ScoreMatrix (int match, int mismatch, int gap);

        /**
         * getScore()
         *
         * @brief
         * Calculate the score for two characters from the haystack and needle.
         * @param r1 Character from haystack.
         * @param r2 Character from needle.
         */
        int
        getScore (char r1, char r2);

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
        FuzzyQuery (Tp haystackSequence, Tp needleSequence);

        /**
         * setItemParser()
         *
         * @brief
         * Set parser for a single item from a match.
         * @param parser The function used for parsing part of a sequence with start and end.
         */
        void
        setItemParser (std::string (*parser) (Item<Tp>&));

        /**
         * initializeScoreMatrix()
         *
         * @brief
         * Set rewards and penalties for score calculation as well as amounts.
         * @param scoreSet The values used to calculate scores for sequence parts.
         * @param amount The amount of matches to be retrieved.
         */
        void
        initializeScoreMatrix (ScoreMatrix scoreSet, int amount);

        /**
         * initializeMatrix()
         *
         * @brief
         * Span the matrix of haystack*needle dimensions and set undefined states.
         * @param haystackSequence The complete genome to be searched.
         * @param needleSequence The needle sequence to be approximately searched.
         */
        void
        initializeMatrix (Tp haystackSequence, Tp needleSequence);

        /**
         * updateMatrix()
         *
         * @brief
         * Update all nodes in the matrix.
         */
        void
        updateMatrix ();

        /**
         * search()
         *
         * @brief
         * Get aligned haystack sequence and needle sequence matches according to amount.
         */
        std::list<Match<Tp>>&
        search ();


        Tp _needleSequence;
        Tp _haystackSequence;
        int _score;
        int _amount;
        ScoreMatrix _scoreSet;
        std::vector<std::vector<Node>> _nodes;
        std::list<Match<Tp>> _matches;
        std::string (*_parser)(Item<Tp>&);
    };

    ScoreMatrix continuityMatrix { 1, 0, 2 }, disparityMatrix { 1, -1, 0 }, mixedMatrix { 1, -1, 1 };
}

typedef sqn::SimpleType<unsigned char, sqn::Dna5Specification> Dna5;
typedef sqn::Sequence<Dna5> Dna5Sequence;

#endif
