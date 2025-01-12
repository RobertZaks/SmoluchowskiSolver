#include "indexlist.hpp"

IndexList::IndexList(unsigned allocsize)
    : alloc_size(0), current_size(0), arr(0), 
    first_not_listed_index(0)
{
    Resize(allocsize);
}

IndexList::~IndexList()
{
    if(arr) {
        delete [] arr;
    }
}

unsigned IndexList::GetFirstNotListed() const
{
    return first_not_listed_index;
}

unsigned IndexList::GetCurrentSize() const
{
    return current_size;
}

unsigned IndexList::operator[](unsigned i) const
{
#if HANDLEXCEPTION
    if(i >= current_size) {
        throw "assume i < current_size";
    }
#endif
    return arr[i];
}

void IndexList::Resize(unsigned need_size)
{
    unsigned *tmp;
    unsigned i;
    if(need_size > alloc_size) {
        alloc_size = need_size;
        tmp = new unsigned[alloc_size];
        for(i = 0; i < current_size; ++i) {
            tmp[i] = arr[i];
        }
        delete [] arr;
        arr = tmp;
    }
}

void IndexList::AddIndex(unsigned index)
{
    unsigned i, j;
    i = first_not_listed_index;
    for(; i < current_size; ++i) {
        if(arr[i] > index) {
            break;
        }
    }
    if(current_size == alloc_size) {
        Resize(alloc_size + indexlist_add_allocsize);
    }
    for(j = current_size; j > i; --j) {
        arr[j] = arr[j - 1];
    }
    arr[i] = index;
    current_size++;
    if(index != first_not_listed_index) {
        return;
    }
    for(j = i; j + 1 < current_size; ++j) {
        if((arr[j] + 1) != arr[j + 1]) {
            break;
        }
    }
    first_not_listed_index = arr[j] + 1;
}


/*
void IndexList::AddIndex(unsigned index)
{
    Item *tmp, *pitem;
    tmp = new Item;
    tmp->index = index;
    if(!first || (first->index > index)) {
        tmp->next = first;
        first = tmp;
    } else {
        pitem = first;
        while(pitem->next && (pitem->next->index < index)) {
            pitem = pitem->next;
        }
        tmp->next = pitem->next;
        pitem->next = tmp;
    }
    if(index != first_not_listed_index) {
        return;
    }
    while(tmp->next && (tmp->index + 1 == tmp->next->index)) 
    {
        tmp = tmp->next;
    }
    first_not_listed_index = tmp->index + 1;
}
*/
