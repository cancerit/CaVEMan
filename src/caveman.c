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

#include <setup.h>
#include <split.h>
#include <mstep.h>
#include <merge.h>
#include <estep.h>
#include <dbg.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

void print_usage (int exit_code){
	printf("Usage: caveman [setup|split|mstep|merge|estep] [-v] [-h] \n\n");
	printf("-v  --version         Print version information\n");
	printf("-h  --help            Display this usage information.\n");
  	exit(exit_code);
}

void setup_options(int argc, char *argv[]){
	const struct option long_opts[] =
	{
             	/*{"setup", no_argument, 0, 's'},
             	{"split", no_argument, 0, 'p'},
             	{"mstep", no_argument, 0, 'm'},
             	{"merge", no_argument, 0, 'e'},
             	{"estep", no_argument, 0, 't'},    */
             	{"help", no_argument, 0, 'h'}, 
             	{"version", no_argument, 0, 'v'},               	        	
             	{ NULL, 0, NULL, 0}
   }; //End of declaring opts
   
   int index = 0;
   int iarg = 0;   
   
   if(argc < 2){
   	print_usage(0);
   }else if(argc >= 2){
   	char *option = argv[1];
		if(strcmp(option,"setup")==0){
			argv[1]= "";
			setup_main(argc,argv);
			return;
		}else if(strcmp(option,"split")==0){
			argv[1]= "";
			split_main(argc,argv);
			return;
		}else if(strcmp(option,"mstep")==0){
			argv[1]= "";
			mstep_main(argc,argv);
			return;
		}else if(strcmp(option,"merge")==0){
			argv[1]= "";
			merge_main(argc,argv);
			return;
		}else if(strcmp(option,"estep")==0){
			argv[1]= "";
			estep_main(argc,argv);
			return;
		}		
   }
         
   //Iterate through options
   while((iarg = getopt_long(argc, argv, "vh", long_opts, &index)) != -1){
   	switch(iarg){
   		case 'h':
   			print_usage(0);
         	break;
   		
			case 'v':
				printf ("%s\n",CAVEMAN_VERSION);
				return;
	
			case '?':
				print_usage(1);
				break;
				
      	default:
      		print_usage(1);
                           
   	}; // End of args switch statement
   	
   }//End of iteration through options   
   return;
}

int main(int argc, char *argv[]){
	setup_options(argc,argv);	  
   return 0;
}