/**   LICENSE
* Copyright (c) 2014 Genome Research Ltd. 
* 
* Author: Cancer Genome Project cgpit@sanger.ac.uk 
* 
* This file is part of CaVEMan. 
* 
* CaVEMan is free software: you can redistribute it and/or modify it under 
* the terms of the GNU Affero General Public License as published by the Free 
* Software Foundation; either version 3 of the License, or (at your option) any 
* later version. 
* 
* This program is distributed in the hope that it will be useful, but WITHOUT 
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
* FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more 
* details. 
* 
* You should have received a copy of the GNU Affero General Public License 
* along with this program. If not, see <http://www.gnu.org/licenses/>. 
*/

#include <stdlib.h>

#ifdef ELEMENT_TYPE
#ifdef ELEMENTS_PER_NODE

#include <templates.h>

#define ListNode_TYPE TEMPLATE(ELEMENT_TYPE,ListNode)
#define List_TYPE TEMPLATE(ELEMENT_TYPE,List)

struct ListNode_TYPE;

typedef struct ListNode_TYPE {
	struct ListNode_TYPE *next;
	struct ListNode_TYPE *prev;
	int count;
	ELEMENT_TYPE values[ELEMENTS_PER_NODE];
} ListNode_TYPE;

typedef struct List_TYPE {
	int count;
	ListNode_TYPE *first;
	ListNode_TYPE *last;
} List_TYPE;

List_TYPE *TEMPLATE(ELEMENT_TYPE,List_create)();
void TEMPLATE(ELEMENT_TYPE,List_destroy)(List_TYPE *list);

#define List_count(A) ((A)->count)
#define List_first(A) ((A)->first->values[0])
#define List_last(A) ((A)->last->values[(A)->last->count-1])
					
void TEMPLATE(ELEMENT_TYPE,List_push)(List_TYPE *list, ELEMENT_TYPE value);
int TEMPLATE(ELEMENT_TYPE,List_pop)(List_TYPE *list, ELEMENT_TYPE *result);

int TEMPLATE(ELEMENT_TYPE,List_unshift)(List_TYPE *list, ELEMENT_TYPE *result);
void TEMPLATE(ELEMENT_TYPE,List_shift)(List_TYPE *list, ELEMENT_TYPE value);

List_TYPE *TEMPLATE(ELEMENT_TYPE,List_copy)(List_TYPE *list);
List_TYPE *TEMPLATE(ELEMENT_TYPE,List_join)(List_TYPE *list1, List_TYPE *list2);
void TEMPLATE(ELEMENT_TYPE,List_split)(List_TYPE *list, int split_index, List_TYPE *left, List_TYPE *right);

void TEMPLATE(ELEMENT_TYPE,List_splitNode)(List_TYPE *list, ListNode_TYPE *node);
void TEMPLATE(ELEMENT_TYPE,List_insert)(List_TYPE *list, ListNode_TYPE *node, int index, ELEMENT_TYPE value);

#define LIST_FOR_EACH_NODE(ELEMENT_TYPE, L, S, M, V) TEMPLATE(ELEMENT_TYPE,ListNode) *_node = NULL; \
	TEMPLATE(ELEMENT_TYPE,ListNode) *V = NULL;			\
	for (V = _node = L->S; _node != NULL; V = _node = _node->M)

#define LIST_FOR_EACH_ELEMENT(ELEMENT_TYPE, L, S, M, E) TEMPLATE(ELEMENT_TYPE,ListNode) *_node = NULL; \
	int _i;								\
	ELEMENT_TYPE E;							\
	for (_node = L->S; _node != NULL; _node = _node->M)		\
		for (E = _node->values[_i = 0]; _i < _node->count; E = _node->values[++_i])

#define LIST_FOR_EACH_ELEMENT_MORE(ELEMENT_TYPE, L, S, M, E, F) TEMPLATE(ELEMENT_TYPE,ListNode) *_node = NULL; \
	int _i;								\
	ELEMENT_TYPE E;							\
	int F;								\
	for (_node = L->S; _node != NULL; _node = _node->M)		\
		for (E = _node->values[_i = 0], F = (_i < _node->count-1) || (_node->next != NULL); \
		     _i < _node->count;					\
		     E = _node->values[++_i], F = (_i < _node->count-1) || (_node->next != NULL))

#define List_clear(ELEMENT_TYPE, A) {LIST_FOR_EACH_ELEMENT(ELEMENT_TYPE, A, first, next, cur) {free(cur);}}
#define List_clear_destroy(ELEMENT_TYPE, A) {List_clear(ELEMENT_TYPE, A); TEMPLATE(ELEMENT_TYPE,List_destroy)(A);}

#undef ListNode_TYPE
#undef List_TYPE

#endif
#endif
