#ifndef INDEXLIST_HPP
#define INDEXLIST_HPP

enum {indexlist_add_allocsize = 3};

/* list of increasing indexes */
/* don't add one index many times! */
class IndexList {
    unsigned alloc_size;
    unsigned current_size;
    unsigned *arr;
    unsigned first_not_listed_index;
public:
    IndexList(unsigned allocsize = 0);
    unsigned GetFirstNotListed() const;
    void AddIndex(unsigned index);
    void Resize(unsigned need_size);
    unsigned GetCurrentSize() const;
    unsigned operator[](unsigned i) const;
    ~IndexList();
};


#endif /* INDEXLIST_HPP */
