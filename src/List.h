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

#ifndef _List_h
#define _List_h

#include <stdlib.h>

#define ELEMENTS_PER_NODE 8

struct ListNode;

typedef struct ListNode{
	struct ListNode *next;
	struct ListNode *prev;
        int numElements;
	void *values[ELEMENTS_PER_NODE];
} ListNode;

typedef struct List{
		int count;
		ListNode *first;
		ListNode *last;
} List;

List *List_create();
void List_destroy(List *list);
void List_clear(List *list);
void List_clear_destroy(List *list);

#define List_count(A) ((A)->count)
#define List_first(A) ((A)->first != NULL ? (A)->first->values[0] : NULL)
#define List_last(A) ((A)->last != NULL ? (A)->last->values[(A)->last->numElements-1] : NULL)

void List_push(List *list, void *value);
void *List_pop(List *list);

void *List_unshift(List *list);
void List_shift(List *list, void *value);

List *List_copy(List *list);
List *List_join(List *list1, List *list2);
void List_split(List *list, int split_index, List *left, List *right);

void List_maybeSplitNode(List *list, ListNode *node);
void List_insert(List *list, ListNode *node, int index, void *value);

#define LIST_FOR_EACH_NODE(L, S, M, V) ListNode *_node = NULL; \
  ListNode *V = NULL;					       \
  for (V = _node = L->S; _node != NULL; V = _node = _node->M)

#define LIST_FOR_EACH_ELEMENT(L, S, M, E) ListNode *_node = NULL;	\
  int _i;								\
  void *E;								\
  for (_node = L->S; _node != NULL; _node = _node->M)			\
    for (E = _node->values[_i = 0]; _i < _node->numElements; E = _node->values[++_i])

#define LIST_FOR_EACH_ELEMENT_MORE(L, S, M, E, F) ListNode *_node = NULL; \
  int _i;								\
  void *E;								\
  int F;								\
  for (_node = L->S; _node != NULL; _node = _node->M)			\
    for (E = _node->values[_i = 0], F = (_i < _node->numElements-1) || (_node->next != NULL); \
	 _i < _node->numElements;					\
	 E = _node->values[++_i], F = (_i < _node->numElements-1) || (_node->next != NULL)) 

#endif
