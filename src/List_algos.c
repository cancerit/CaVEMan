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

#include <List_algos.h>
#include <dbg.h>

/*
inline void ListNode_swap(ListNode *a, ListNode *b)
{
    void *temp = a->value;
    a->value = b->value;
    b->value = temp;
}
*/

int List_bubble_sort(List *list, List_compare cmp)
{
    int sorted = 1;

    if(List_count(list) <= 1) {
        return 0;  // already sorted
    }

    do {
        sorted = 1;
        LIST_FOREACH(list, first, next, cur) {
	  int curi;
	  for (curi=0; curi<cur->numElements; ++curi) {
	    if (curi<cur->numElements-1) {
	      if (cmp(cur->values[curi], cur->values[curi+1]) > 0) {
		void* tmp = cur->values[curi];
		cur->values[curi] = cur->values[curi+1];
		cur->values[curi+1] = tmp;
		sorted = 0;
	      }
	    } else if(cur->next) {
	      if(cmp(cur->values[curi], cur->next->values[0]) > 0) {
		void* tmp = cur->values[curi];
		cur->values[curi] = cur->next->values[0];
		cur->next->values[0] = tmp;
		sorted = 0;
	      }
            }
	  }
        }
    } while(!sorted);

    return 0;
}

inline List *List_merge(List *left, List *right, List_compare cmp)
{
    List *result = List_create();
    void *val = NULL;

    while(List_count(left) > 0 || List_count(right) > 0) {
        if(List_count(left) > 0 && List_count(right) > 0) {
            if(cmp(List_first(left), List_first(right)) <= 0) {
                val = List_unshift(left);
            } else {
                val = List_unshift(right);
            }

            List_push(result, val);
        } else if(List_count(left) > 0) {
            val = List_unshift(left);
            List_push(result, val);
        } else if(List_count(right) > 0) {
            val = List_unshift(right);
            List_push(result, val);
        }
    }

    return result;
}

List *List_merge_sort(List *list, List_compare cmp)
{
    if(List_count(list) <= 1) {
        return list;
    }

    List *left = List_create();
    List *right = List_create();
    int middle = List_count(list) / 2;

    LIST_FOREACH(list, first, next, cur) {
      int curi;
      for (curi=0; curi<cur->numElements; ++curi) {
        if(middle > 0) {
            List_push(left, cur->values[curi]);
        } else {
            List_push(right, cur->values[curi]);
        }

        middle--;
      }
    }

    List *sort_left = List_merge_sort(left, cmp);
    List *sort_right = List_merge_sort(right, cmp);

    if(sort_left != left) List_destroy(left);
    if(sort_right != right) List_destroy(right);

    return List_merge(sort_left, sort_right, cmp);
}
