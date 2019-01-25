/**   LICENSE
* Copyright (c) 2014-2018 Genome Research Ltd.
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
*
*    1. The usage of a range of years within a copyright statement contained within
*    this distribution should be interpreted as being equivalent to a list of years
*    including the first and last year specified and all consecutive years between
*    them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
*    2009, 2011-2012’ should be interpreted as being identical to a statement that
*    reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
*    statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
*    identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
*    2009, 2010, 2011, 2012’."
*
*/

#include <List.h>
#include <dbg.h>
#include <assert.h>

List *List_create(){
	return calloc(1, sizeof(List));
}

void List_destroy(List *list){
	assert(list != NULL);
	LIST_FOREACH(list, first, next, cur){
		if(cur->prev){
			free(cur->prev);
		}
	}
	free(list->last);
	free(list);
}

void List_clear(List *list){
	assert(list != NULL);
	LIST_FOREACH(list, first, next, cur){
		if(cur->value) free(cur->value);
	}
}

void List_clear_destroy(List *list){
	assert(list != NULL);
	if(List_count(list) > 0){
		List_clear(list);
	}
	List_destroy(list);
}

void List_push(List *list, void *value){
	assert(list != NULL);
	ListNode *node = calloc(1, sizeof(ListNode));
	check_mem(node);

	node->value = value;

	if(list->last == NULL){
		list->first = node;
		list->last = node;
	}else{
		list->last->next = node;
		node->prev = list->last;
		list->last = node;
	}

	list->count++;

error:
	return;
}

void *List_pop(List *list){
	assert(list != NULL);
	ListNode *node = list->last;
	return node != NULL ? List_remove(list, node) : NULL;
}

void List_shift(List *list, void *value)
{
    ListNode *node = calloc(1, sizeof(ListNode));
    check_mem(node);

    node->value = value;

    if(list->first == NULL) {
        list->first = node;
        list->last = node;
    } else {
        node->next = list->first;
        list->first->prev = node;
        list->first = node;
    }

    list->count++;

error:
    return;
}

void *List_unshift(List *list)
{
    ListNode *node = list->first;
    return node != NULL ? List_remove(list, node) : NULL;
}

List *List_copy(List *list){
	assert(list != NULL);
	List *listCopy = List_create();
	LIST_FOREACH(list, first, next, cur){
		List_push(listCopy, cur->value);
	}
	return listCopy;
}

List *List_join(List *list1, List *list2){
	assert(list1 != NULL);
	assert(list2 != NULL);
	List *joined = List_copy(list1);
	LIST_FOREACH(list2, first, next, cur){
		List_push(joined, cur->value);
	}
	return joined;
}

void List_split(List *list, int split_pos, List *left, List *right){
	assert(list != NULL);
	assert(List_count(list)>0);
	assert(split_pos <= List_count(list));
	LIST_FOREACH(list, first, next, cur){
		if(split_pos > 0){
			List_push(left,&cur->value);
		}else{
			List_push(right,&cur->value);
		}
		split_pos--;
	}
	return;
}

void *List_remove(List *list, ListNode *node){
	assert(list != NULL);
	void *result = NULL;

	check(list->first && list->last, "List is empty.");
	check(node, "node can't be NULL");

	if(node == list->first && node == list->last){
		list->first = NULL;
		list->last = NULL;
	}else if(node == list->first){
		list->first = node->next;
		check(list-> first != NULL, "Invalid list, somehow got a first that is NULL.");
		list->first->prev = NULL;
	}else if(node == list->last){
		list->last = node->prev;
		check(list->last != NULL, "Invalid list, somehow got a last that is NULL.");
		list->last->next = NULL;
	}else{
		ListNode *after = node->next;
		ListNode *before = node->prev;
		after->prev = before;
		before->next = after;
	}

	list->count--;
	result = node->value;
	free(node);

error:
	return result;
}
