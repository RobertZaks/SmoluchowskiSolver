#ifndef INDEXLIST_HPP
#define INDEXLIST_HPP

/*!
    \file indexlist.hpp
    \brief List of increasing indexes

    The file defines the class of list of increasing indexes which allow
    you to get first not listed index.
*/



//! Number of size which added to indexlist if needed
enum {indexlist_add_allocsize = 3};

//! Class of increasing indexes list
class IndexList {
    unsigned alloc_size;
    unsigned current_size;
    unsigned *arr;
    unsigned first_not_listed_index;
    // Add size to list if needeed
    void Resize(unsigned need_size);
public:
    //! Initialize list with allocsize
    IndexList(unsigned allocsize = 0);
    //! Get first index wich not listed in list
    unsigned GetFirstNotListed() const;
    //! Add index to list; Be carefull: don't add one index many times!
    void AddIndex(unsigned index);
    //! Get number of saved indeces
    unsigned GetCurrentSize() const;
    //! Get i'th index in increasing order
    unsigned operator[](unsigned i) const;
    //! The destructor
    ~IndexList();
};


#endif /* INDEXLIST_HPP */
