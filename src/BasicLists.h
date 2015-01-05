/**   LICENSE
* Copyright (c) 2015 Genome Research Ltd. 
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
