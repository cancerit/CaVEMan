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

#include <assert.h>

#define List_TYPE TEMPLATE(ELEMENT_TYPE,List)
#define List_COMPARE TEMPLATE(ELEMENT_TYPE,List_compare)

typedef int (*List_COMPARE)(const ELEMENT_TYPE a, const ELEMENT_TYPE b);

int TEMPLATE(ELEMENT_TYPE,List_bubble_sort)(List_TYPE *list, List_COMPARE cmp);

List_TYPE *TEMPLATE(ELEMENT_TYPE,List_merge_sort)(List_TYPE *list, List_COMPARE cmp);

List_TYPE *TEMPLATE(ELEMENT_TYPE,List_merge)(List_TYPE *left, List_TYPE *right, List_COMPARE cmp);

#undef List_TYPE
#undef List_COMPARE

#endif
