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

#include <List.h>
#include <dbg.h>
#include <assert.h>

List *List_create() {
  return calloc(1, sizeof(List));
}

void List_destroy(List *list) {
  assert(list != NULL);
  LIST_FOR_EACH_NODE(list, first, next, cur) {
    if(cur->prev){
      free(cur->prev);
    }
  }
  free(list->last);
  free(list);
}

void List_clear(List *list) {
  assert(list != NULL);
  LIST_FOR_EACH_ELEMENT(list, first, next, cur) {
    free(cur);
  }
}

void List_clear_destroy(List *list) {
  assert(list != NULL);
  if(List_count(list) > 0){
    List_clear(list);
  }	
  List_destroy(list);
}

void List_push(List *list, void *value) {
  assert(list != NULL);
  if ((list->last == NULL) || (list->last->numElements == ELEMENTS_PER_NODE)) {
    ListNode *node = calloc(1, sizeof(ListNode));
    check_mem(node);
    node->numElements = 1;
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
    list->last->values[list->last->numElements++] = value;
  }
  list->count++;
  
error:
  return;
}

void *List_pop(List *list) {
  assert(list != NULL);
  ListNode *node = list->last;
  if (node != NULL) {
    void *result = node->values[--node->numElements];
    if (node->numElements == 0) {
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
    return result;
  } else {
    return NULL;
  }
}

void List_shift(List *list, void *value) {
  assert(list != NULL);
  if ((list->first == NULL) || (list->first->numElements == ELEMENTS_PER_NODE)) {
    ListNode *node = calloc(1, sizeof(ListNode));
    check_mem(node);
    node->numElements = 1;
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
    for (i=list->first->numElements++; i>0; --i) {
      list->first->values[i] = list->first->values[i-1];
    }
    list->first->values[0] = value;
  }
  list->count++;

error:
  return;
}

void *List_unshift(List *list) {
  assert(list != NULL);
  ListNode *node = list->first;
  if (node != NULL) {
    void *result = node->values[0];
    if (--node->numElements == 0) {
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
      for (i=0; i<node->numElements; ++i) {
	node->values[i] = node->values[i+1];
      }
    }
    list->count--;
    return result;
  } else {
    return NULL;
  }
}

List *List_copy(List *list) {
  assert(list != NULL);
  List *listCopy = List_create();
  LIST_FOR_EACH_ELEMENT(list, first, next, cur) {
      List_push(listCopy, cur);
  }
  return listCopy;
}

List *List_join(List *list1, List *list2) {
  assert(list1 != NULL);
  assert(list2 != NULL);
  List *joined = List_copy(list1);
  LIST_FOR_EACH_ELEMENT(list2, first, next, cur) {
      List_push(joined, cur);
  }
  return joined;
}

void List_split(List *list, int split_pos, List *left, List *right) {
  assert(list != NULL);
  assert(List_count(list)>0);
  assert(split_pos <= List_count(list));
  LIST_FOR_EACH_ELEMENT(list, first, next, cur) {
      if(split_pos > 0) {
	List_push(left,&cur);
      } else {
	List_push(right,&cur);
      }
      split_pos--;
  }
  return;
}

void List_maybeSplitNode(List *list, ListNode *node) {
  if (node->numElements == ELEMENTS_PER_NODE) {
    ListNode *newNode = calloc(1, sizeof(ListNode));
    newNode->next = node->next;
    newNode->prev = node;
    newNode->numElements = ELEMENTS_PER_NODE/2;
    int i;
    for (i=0; i<ELEMENTS_PER_NODE/2; ++i) {
      newNode->values[i] = node->values[i+ELEMENTS_PER_NODE/2];
    }
    node->next = newNode;
    node->numElements -= ELEMENTS_PER_NODE/2;
    if (list->last == node) list->last = newNode;
  }
}

void List_insert(List *list, ListNode *node, int index, void *value) {
  List_maybeSplitNode(list, node);
  if (index < node->numElements) {
    int i;
    for (i=++node->numElements; i>index; --i) {
      node->values[i] = node->values[i-1];
    }
    node->values[index] = value;
    list->count++;
  } else {
    List_insert(list, node->next, index-node->numElements, value);
  }
}
