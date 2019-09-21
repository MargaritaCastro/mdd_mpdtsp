/*
 * -----------------------------------------------------------
 * Reversible linked list
 * -----------------------------------------------------------
 */

#ifndef LINKEDLIST_HPP_
#define LINKEDLIST_HPP_

#include <cstdlib>
#include <ilcp/cpext.h>

/**
 * Element of a reversible linked list
 */
struct LinkedListElem {

    IlcRevAny   next;       /**< next element in list */

    /**
     * Constructor
     */
    LinkedListElem(IloCP cp) {
        next.setValue(cp, NULL);
    }
};


/**
 * Reversible Doubly Linked List
 */
template<class ListType>
struct RevLinkedList {

    IloCP                   cp;             /**< solver object */
    ListType                *header;    /**< first dummy element in the list */
    int                     max_id;         /**< maximum acceptable id */

    /**
     * Empty constructor
     */
    RevLinkedList() { }

    /**
     * Initialize list.
     * Requires context to define first element
     **/
    void initialize(IloCP _cp, ListType* _header, int _max_id) {
        assert( _header != NULL );
        assert( _max_id >= 0 );
        cp =_cp;
        header = _header;
        max_id = _max_id;
        header->next.setValue(cp, NULL);
    }

    /**
     * Clear list
     */
    void clear() {
        header->next.setValue(cp, NULL);
    }

    /**
     * Get first element
     */
    ListType* get_first() {
        return( static_cast<ListType*>(header->next.getValue()) );
    }

    /**
     * Get next element
     */
    ListType* get_next(ListType *type) {
        return( static_cast<ListType*>(type->next.getValue()) );
    }

    /**
     * Extend list with another element. 'Previous' should ideally be the last element
     * of the list. Returns the 'next' element to facilitate iterating on consecutive
     * additions
     */
    ListType* extend(ListType *previous, ListType *next) {
        previous->next.setValue(cp, next);
        return next;
    }

    /**
     * Close list at element, i.e. element is the last of the list
     */
    void close(ListType* elem) {
        elem->next.setValue(cp, NULL);
    }

    /**
     * Get header of the list
     */
    ListType* get_header() {
        return header;
    }

    /**
     * Check if list is empty
     */
    bool is_empty() {
        return( header->next.getValue() == NULL );
    }

    /**
     * Check if list has a single element
     */
    bool is_unique() {
    	if( header->next.getValue() == NULL )
    		return false;
    	return( (static_cast<ListType*>(header->next.getValue()))->next.getValue() == NULL );
    }
};



#endif /* LINKEDLIST_HPP_ */
