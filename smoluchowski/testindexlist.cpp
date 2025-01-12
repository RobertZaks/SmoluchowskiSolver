#include "indexlist.hpp"
#include <stdio.h>

void print_indexlist(const IndexList &list)
{
    unsigned i, size;
    size = list.GetCurrentSize();
    for(i = 0; i < size; ++i) {
        printf("%d ", list[i]);
    }
    putchar('\n');
}

int main()
{
    IndexList index_ignore;
    printf("first_not_listed: %d\n", index_ignore.GetFirstNotListed());
    index_ignore.AddIndex(5);
    print_indexlist(index_ignore);
    printf("first_not_listed: %d\n", index_ignore.GetFirstNotListed());
    index_ignore.AddIndex(1);
    print_indexlist(index_ignore);
    printf("first_not_listed: %d\n", index_ignore.GetFirstNotListed());
    index_ignore.AddIndex(3);
    print_indexlist(index_ignore);
    printf("first_not_listed: %d\n", index_ignore.GetFirstNotListed());
    index_ignore.AddIndex(4);
    print_indexlist(index_ignore);
    printf("first_not_listed: %d\n", index_ignore.GetFirstNotListed());
    index_ignore.AddIndex(0);
    print_indexlist(index_ignore);
    printf("first_not_listed: %d\n", index_ignore.GetFirstNotListed());
    index_ignore.AddIndex(6);
    print_indexlist(index_ignore);
    printf("first_not_listed: %d\n", index_ignore.GetFirstNotListed());
    /* don't do this because we want to have first listed index in 
     * index_ignore as element in index_ignore on first not listed index
    index_ignore.AddIndex(2);
    print_indexlist(index_ignore);
    printf("first_not_listed: %d\n", index_ignore.GetFirstNotListed());
    */
    index_ignore.AddIndex(2);
    print_indexlist(index_ignore);
    printf("first_not_listed: %d\n", index_ignore.GetFirstNotListed());
    return 0;
}
