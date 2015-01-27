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

#include <templates.h>

#include <List_algos.h>
#include <dbg.h>

#define ListNode_TYPE TEMPLATE(ELEMENT_TYPE,ListNode)
#define List_TYPE TEMPLATE(ELEMENT_TYPE,List)
#define List_COMPARE TEMPLATE(ELEMENT_TYPE,List_compare)

int TEMPLATE(ELEMENT_TYPE,List_bubble_sort)(List_TYPE *list, List_COMPARE cmp) {
	int sorted = 1;

	if(List_count(list) <= 1) {
		return 0;  // already sorted
	}

	do {
		sorted = 1;
		LIST_FOR_EACH_NODE(ELEMENT_TYPE, list, first, next, cur) {
			int curi;
			for (curi=0; curi<cur->count; ++curi) {
				if (curi<cur->count-1) {
					if (cmp(cur->values[curi], cur->values[curi+1]) > 0) {
						ELEMENT_TYPE tmp = cur->values[curi];
						cur->values[curi] = cur->values[curi+1];
						cur->values[curi+1] = tmp;
						sorted = 0;
					}
				} else if(cur->next) {
					if(cmp(cur->values[curi], cur->next->values[0]) > 0) {
						ELEMENT_TYPE tmp = cur->values[curi];
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

inline List_TYPE *TEMPLATE(ELEMENT_TYPE,List_merge)(List_TYPE *left, List_TYPE *right, List_COMPARE cmp) {
	List_TYPE *result = TEMPLATE(ELEMENT_TYPE,List_create)();
	ELEMENT_TYPE val;

	while(List_count(left) > 0 || List_count(right) > 0) {
		if(List_count(left) > 0 && List_count(right) > 0) {
			if(cmp(List_first(left), List_first(right)) <= 0) {
				TEMPLATE(ELEMENT_TYPE,List_unshift)(left, &val);
			} else {
				TEMPLATE(ELEMENT_TYPE,List_unshift)(right, &val);
			}

			TEMPLATE(ELEMENT_TYPE,List_push)(result, val);
		} else if(List_count(left) > 0) {
			TEMPLATE(ELEMENT_TYPE,List_unshift)(left, &val);
			TEMPLATE(ELEMENT_TYPE,List_push)(result, val);
		} else if(List_count(right) > 0) {
			TEMPLATE(ELEMENT_TYPE,List_unshift)(right, &val);
			TEMPLATE(ELEMENT_TYPE,List_push)(result, val);
		}
	}

	return result;
}

List_TYPE *TEMPLATE(ELEMENT_TYPE,List_merge_sort)(List_TYPE *list, List_COMPARE cmp) {
	if(List_count(list) <= 1) {
		return list;
	}

	List_TYPE *left = TEMPLATE(ELEMENT_TYPE,List_create)();
	List_TYPE *right = TEMPLATE(ELEMENT_TYPE,List_create)();
	int middle = List_count(list) / 2;

	LIST_FOR_EACH_ELEMENT(ELEMENT_TYPE, list, first, next, cur) {
		if(middle > 0) {
			TEMPLATE(ELEMENT_TYPE,List_push)(left, cur);
		} else {
			TEMPLATE(ELEMENT_TYPE,List_push)(right, cur);
		}

		middle--;
	}

	List_TYPE *sort_left = TEMPLATE(ELEMENT_TYPE,List_merge_sort)(left, cmp);
	List_TYPE *sort_right = TEMPLATE(ELEMENT_TYPE,List_merge_sort)(right, cmp);

	if(sort_left != left) TEMPLATE(ELEMENT_TYPE,List_destroy)(left);
	if(sort_right != right) TEMPLATE(ELEMENT_TYPE,List_destroy)(right);

	return TEMPLATE(ELEMENT_TYPE,List_merge)(sort_left, sort_right, cmp);
}

#undef ListNode_TYPE
#undef List_TYPE
#undef List_COMPARE

#endif
