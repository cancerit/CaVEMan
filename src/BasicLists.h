#ifndef _BasicLists_h
#define _BasicLists_h

#define ELEMENT_TYPE int
#define ELEMENTS_PER_NODE 16
#include <List.h>
#undef ELEMENT_TYPE
#undef ELEMENTS_PER_NODE

#define ELEMENT_TYPE float
#define ELEMENTS_PER_NODE 16
#include <List.h>
#undef ELEMENT_TYPE
#undef ELEMENTS_PER_NODE

typedef char* String;

#define ELEMENT_TYPE String
#define ELEMENTS_PER_NODE 8
#include <List.h>
#undef ELEMENT_TYPE
#undef ELEMENTS_PER_NODE

#endif
