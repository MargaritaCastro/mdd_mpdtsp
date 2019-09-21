/*
 * -----------------------------------------------------------
 * Data structure to store positive sets of integers.

 * Should be used in the cases where one needs to do
 * fast intersections and unions, and does not iterate
 * on all elements too many times.
 * -----------------------------------------------------------
 */

#ifndef INTSET_HPP_
#define INTSET_HPP_

#define NOT_COMPUTED -1     /**< indicates if the size was not computed */

#include <boost/dynamic_bitset.hpp>
#include <cassert>
#include <iostream>
#include <fstream>
#include <ilcp/cpext.h>
#include "../util/util.hpp"


/**
 * Integer Set structure
 */
struct IntSet {

    /** Constructor */
    IntSet(int _min, int _max, bool _filled);

    /** Empty constructor */
    IntSet();

    /** Check if set contains element */
    bool contains(int elem);

    /** Add an element to the set */
    void add(int elem);

    /** Add all possible elements to the set */
    void add_all_elements();

    /** Remove element, if it is contained */
    void remove(int elem);

    /** Get number of elements in the set */
    int get_size();

    /** Get the first element of the set */
    int get_first();

    /** Get next element higher than the one passed as parameter */
    int get_next(int elem);

    /** Get end of the set (beyond last element) */
    const int get_end();

    /** Clear set */
    void clear();

    /** Resize */
    void resize(int _min, int _max, bool _filled);

    /** Add all elements between min and max */
    void fill();

    /** Take the union with another intset */
    void union_with(IntSet& intset);

    /** Take the intersection with another intset */
    void intersect_with(IntSet& intset);

    /** Set minus operation */
    void set_minus(IntSet& rhs);

    /** Assignment operator */
    IntSet& operator=(const IntSet& rhs);

    /** Comparison operators */
    bool operator == (const IntSet& intset);
    bool operator < (const IntSet& intset);


    // parameters

    boost::dynamic_bitset<>     set;            /**< bitvector representing the set */
    const int                   end;            /**< position beyond end of the set */
    int                         size;           /**< number of elements in the set */
    int                         min;            /**< minimum possible element of the set */
    int                         max;            /**< maximum possible element of the set */
//    int                         shift;          /**< shift of element to be added in the set */
};


/**
 * Lexicographic comparator function for IntSet class.
 */
struct IntSetLexLessThan {
    bool operator()(const IntSet* setA, const IntSet* setB) const {
        return setA->set < setB->set;
    }
};



/**
 * -----------------------------------------------
 * Inline implementations
 * -----------------------------------------------
 */

/**
 * Constructor
 */
inline IntSet::IntSet(int _min, int _max, bool _filled) : end((int)set.npos) {
    resize(_min, _max, _filled);
    size = NOT_COMPUTED;
}

/**
 * Empty constructor
 */
inline IntSet::IntSet() : end((int)set.npos) {
}


/**
 * Add an element to the set
 */
inline bool IntSet::contains(int elem) {
    assert( elem >= min && elem <= max );
    return( set.test(elem) );
}


/**
 * Add an element to the set
 */
inline void IntSet::add(int elem) {
    assert( elem >= min && elem <= max );
    set.set(elem, true);
    //size = NOT_COMPUTED;
}

/** Remove element, if it is contained */
inline void IntSet::remove(int elem) {
    assert( elem >= min && elem <= max );
    set.set(elem, false);
    //size = NOT_COMPUTED;
}

/**
 * Get the first element of the set
 */
inline int IntSet::get_first() {
    return (set.find_first());
}

/**
 * Get next element higher than the one passed as parameter
 */
inline int IntSet::get_next(int elem) {
    assert( elem >= min && elem <= max );
    for( int v = elem+1; v <= max; v++ ) {
        if( set.test(v) )
            return v;
    }
    return end;
    //return (set.find_next(elem));
}

/**
 * Get end of the set (beyond last element)
 */
inline const int IntSet::get_end() {
    return end;
}

/**
 * Clear set
 */
inline void IntSet::clear() {
    set.reset();
    //size = 0;
}


/**
 * Assignment operator
 */
inline IntSet& IntSet::operator=(const IntSet& rhs) {
    assert(rhs.max == max && rhs.min == min);
    if (this != &rhs) {
        set = rhs.set;
        //size = NOT_COMPUTED;
    }
    return *this;
}


/**
 * Resize
 */
inline void IntSet::resize(int _min, int _max, bool _filled) {

    assert( _min == 0 );

    min = _min;
    max = _max;

    //shift = (-1) * _min;
    set.resize(max - min + 1);

    if( _filled )
        set.set();
    else
        set.reset();

    //size = NOT_COMPUTED;
}

/**
 * Take the union with another intset
 */
inline void IntSet::union_with(IntSet& intset) {
    set |= intset.set;
    //size = NOT_COMPUTED;
}

/**
 * Take the intersection with another intset
 */
inline void IntSet::intersect_with(IntSet& intset) {
    set &= intset.set;
    //size = NOT_COMPUTED;
}


/**
 * Set minus operation
 */
inline void IntSet::set_minus(IntSet& rhs) {
    set -= rhs.set;
    //size = NOT_COMPUTED;
}


/** Get number of elements in the set */
inline int IntSet::get_size() {
//    if( size == NOT_COMPUTED ) {
//        size = set.count();
//    }
//    return size;
    return( set.count() );
}

/**
 * Add all possible elements to the set
 */
inline void IntSet::add_all_elements() {
    set.set();
    //size = max - min + 1;
}

/**
 * Stream output function
 */
inline std::ostream& operator<<(std::ostream &os, IntSet &intset) {
    os << "[ ";
    int val = intset.get_first();
    while( val != intset.get_end() ) {
        os << val << " ";
        val = intset.get_next(val);
    }
    os << "]";
    return os;
}

/**
 * Comparison operators
 */
inline bool IntSet::operator == (const IntSet& intset) {
    return( this->set == intset.set );
}

inline bool IntSet::operator < (const IntSet& intset) {
    return( this->set < intset.set );
}

inline bool operator<(const IntSet& a, const IntSet& b) {
    return( a.set < b.set );
}

inline bool operator==(const IntSet& a, const IntSet& b) {
    return( a.set == b.set );
}



/**
 * Element of a intset linked list
 */
struct IntSetElement {
    int             elem;       /**< integer element */
    IlcRevAny       prev;       /**< previous element */
    IlcRevAny       next;       /**< next element */

    /** Constructor */
    IntSetElement(IloCP cp, int _elem, void* _prev, void* _next)
        : elem(_elem)
    {
        prev.setValue(cp, _prev);
        next.setValue(cp, _next);
    }
};

/**
 * Integer set as a linked list
 */
struct IntSetLinkedList {

    IlcRevAny                first;         /**< first list element */
    IlcRevAny                last;          /**< last list element */
    IloCP                    cp;            /**< solver object */
    int                      max;           /**< maximum number of elements in the list */


    vector<IntSetElement*>   buffer;        /**< buffer of available intset elements */
    IlcRevInt                size;          /**< number of elements */

    /** Constructor */
    IntSetLinkedList(IloCP _cp, int max);

    /** Add element to the list */
    IntSetElement* add(int elem);

    /** Clear list */
    void clear();

    /** Delete element */
    void remove(IntSetElement* elem);
};

/**
 * -----------------------------------------------
 * Inline implementations
 * -----------------------------------------------
 */

/**
 * Constructor of linked list
 */
inline IntSetLinkedList::IntSetLinkedList(IloCP _cp, int _max) : cp(_cp), max(_max) {
    // create header element
    IntSetElement* elem = new( cp.getHeap() ) IntSetElement(cp, -1, NULL, NULL);
    first.setValue(cp, elem);
    last.setValue(cp, elem);

    buffer.resize(max+1);
    for( int i = 0; i < max+1; i++ ) {
        buffer[i] = new( cp.getHeap() ) IntSetElement(cp, i, NULL, NULL);
    }
    size.setValue(cp, 0);
}

/**
 * Add an element to linked list
 */
inline IntSetElement* IntSetLinkedList::add(int _elem) {
    //size.setValue(cp, size.getValue()+1);
    IntSetElement* elem = buffer[_elem];
    //buffer_size.setValue(cp, buffer_size.getValue()-1);

    IntSetElement* last_elem =  static_cast<IntSetElement*>(last.getValue());

    last_elem->next.setValue(cp, elem);
    elem->prev.setValue(cp, last_elem);
    //elem->elem = _elem;
    last.setValue(cp, elem);

    return elem;
}

/**
 * Clear linked list
 */
inline void IntSetLinkedList::clear() {
    last.setValue(cp, first.getValue());
    //buffer_size.setValue(cp, max);
}

/**
 * Delete element
 */
inline void IntSetLinkedList::remove(IntSetElement* elem) {
    size.setValue(cp, size.getValue()-1);
    IntSetElement* prev = ((IntSetElement*)elem->prev.getValue());
    prev->next.setValue(cp, elem->next);
    ((IntSetElement*)elem->next.getValue())->prev.setValue(cp, prev);
}


#endif /* INTSET_HPP_ */

