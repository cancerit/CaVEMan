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

#ifdef ELEMENT_TYPE
#ifdef ELEMENTS_PER_NODE

#include <templates.h>

#include <dbg.h>
#include <assert.h>

#define ListNode_TYPE TEMPLATE(ELEMENT_TYPE,ListNode)
#define List_TYPE TEMPLATE(ELEMENT_TYPE,List)

List_TYPE *TEMPLATE(ELEMENT_TYPE,List_create)() {
	return calloc(1, sizeof(List_TYPE));
}

void TEMPLATE(ELEMENT_TYPE,List_destroy)(List_TYPE *list) {
	assert(list != NULL);
	LIST_FOR_EACH_NODE(ELEMENT_TYPE, list, first, next, cur) {
		if(cur->prev){
			free(cur->prev);
		}
	}
	free(list->last);
	free(list);
}

void TEMPLATE(ELEMENT_TYPE,List_push)(List_TYPE *list, ELEMENT_TYPE value) {
	assert(list != NULL);
	if ((list->last == NULL) || (list->last->count == ELEMENTS_PER_NODE)) {
		ListNode_TYPE *node = calloc(1, sizeof(ListNode_TYPE));
		check_mem(node);
		node->count = 1;
		node->values[0] = value;
		if (list->last == NULL) {
			list->first = node;
			list->last = node;
		} else {
			list->last->next = node;
			node->prev = list->last;
			list->last = node;
		}
	} else {
		list->last->values[list->last->count++] = value;
	}
	list->count++;
  
 error:
	return;
}

int TEMPLATE(ELEMENT_TYPE,List_pop)(List_TYPE *list, ELEMENT_TYPE *result) {
	assert(list != NULL);
	ListNode_TYPE *node = list->last;
	if (node != NULL) {
		*result = node->values[--node->count];
		if (node->count == 0) {
			if (list->first == list->last) {
				list->first = NULL;
				list->last = NULL;
			} else {
				list->last = list->last->prev;
				list->last->next = NULL;
			}
			free(node);
		}
		list->count--;
		return 1;
	} else return 0;
}

void TEMPLATE(ELEMENT_TYPE,List_shift)(List_TYPE *list, ELEMENT_TYPE value) {
	assert(list != NULL);
	if ((list->first == NULL) || (list->first->count == ELEMENTS_PER_NODE)) {
		ListNode_TYPE *node = calloc(1, sizeof(ListNode_TYPE));
		check_mem(node);
		node->count = 1;
		node->values[0] = value;
		if(list->first == NULL) {
			list->first = node;
			list->last = node;
		} else {
			node->next = list->first;
			list->first->prev = node;
			list->first = node;
		}
	} else {
		int i;
		for (i=list->first->count++; i>0; --i) {
			list->first->values[i] = list->first->values[i-1];
		}
		list->first->values[0] = value;
	}
	list->count++;

 error:
	return;
}

int TEMPLATE(ELEMENT_TYPE,List_unshift)(List_TYPE *list, ELEMENT_TYPE *result) {
	assert(list != NULL);
	ListNode_TYPE *node = list->first;
	if (node != NULL) {
		*result = node->values[0];
		if (--node->count == 0) {
			if (list->first == list->last) {
				list->first = NULL;
				list->last = NULL;
			} else {
				list->first = list->first->next;
				list->first->prev = NULL;
			}
			free(node);
		} else {
			int i;
			for (i=0; i<node->count; ++i) {
				node->values[i] = node->values[i+1];
			}
		}
		list->count--;
		return 1;
	} else return 0;
}

List_TYPE *TEMPLATE(ELEMENT_TYPE,List_copy)(List_TYPE *list) {
	assert(list != NULL);
	List_TYPE *listCopy = TEMPLATE(ELEMENT_TYPE,List_create)();
	LIST_FOR_EACH_ELEMENT(ELEMENT_TYPE, list, first, next, cur) {
		TEMPLATE(ELEMENT_TYPE,List_push)(listCopy, cur);
	}
	return listCopy;
}

List_TYPE *TEMPLATE(ELEMENT_TYPE,List_join)(List_TYPE *list1, List_TYPE *list2) {
	assert(list1 != NULL);
	assert(list2 != NULL);
	List_TYPE *joined = TEMPLATE(ELEMENT_TYPE,List_copy)(list1);
	LIST_FOR_EACH_ELEMENT(ELEMENT_TYPE, list2, first, next, cur) {
		TEMPLATE(ELEMENT_TYPE,List_push)(joined, cur);
	}
	return joined;
}

void TEMPLATE(ELEMENT_TYPE,List_split)(List_TYPE *list, int split_pos, List_TYPE *left, List_TYPE *right) {
	assert(list != NULL);
	assert(List_count(list)>0);
	assert(split_pos <= List_count(list));
	LIST_FOR_EACH_ELEMENT(ELEMENT_TYPE, list, first, next, cur) {
		if(split_pos > 0) {
			TEMPLATE(ELEMENT_TYPE,List_push)(left,cur);
		} else {
			TEMPLATE(ELEMENT_TYPE,List_push)(right,cur);
		}
		split_pos--;
	}
	return;
}

void TEMPLATE(ELEMENT_TYPE,List_splitNode)(List_TYPE *list, ListNode_TYPE *node) {
	assert(node->count == ELEMENTS_PER_NODE);
	ListNode_TYPE *newNode = calloc(1, sizeof(ListNode_TYPE));
	if (node->next) {
		newNode->next = node->next;
		newNode->next->prev = newNode;
	} else {
		newNode->next = NULL;
	}
	newNode->prev = node;
	newNode->count = ELEMENTS_PER_NODE/2;
	int i;
	for (i=0; i<ELEMENTS_PER_NODE/2; ++i) {
		newNode->values[i] = node->values[i+ELEMENTS_PER_NODE/2];
	}
	node->next = newNode;
	node->count -= ELEMENTS_PER_NODE/2;
	if (list->last == node) list->last = newNode;
}

void TEMPLATE(ELEMENT_TYPE,List_insert)(List_TYPE *list, ListNode_TYPE *node, int index, ELEMENT_TYPE value) {
	if (index <= node->count) {
		if (node->count == ELEMENTS_PER_NODE)
			TEMPLATE(ELEMENT_TYPE,List_splitNode)(list, node);
		if (index <= node->count) {
			assert(node->count < ELEMENTS_PER_NODE);
			int i;
			for (i=node->count; i>index; --i) {
				node->values[i] = node->values[i-1];
			}
			node->values[index] = value;
			node->count++;
			list->count++;
			return;
		}
	}
	TEMPLATE(ELEMENT_TYPE,List_insert)(list, node->next, index-node->count, value);
}

#undef ListNode_TYPE
#undef List_TYPE

#endif
#endif
