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

#include <covs_access.h>
#include <assert.h>
#include <math.h>
#include <dbg.h>
#include <stdlib.h>

int ********covs_access_generate_cov_array_given_dimensions(int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8){
	int ********array = (int ********)malloc(sizeof(int *******) * dim1);
	check_mem(array);
	int i,j,k,m,n,p,r,s;
	for(i=0;i<dim1;i++){
		array[i] = (int *******)malloc(sizeof(int******) * dim2);
		check_mem(array[i]);
		for(j=0;j<dim2;j++){
			array[i][j] = (int ******)malloc(sizeof(int*****) * dim3);
			check_mem(array[i][j]);
			for(k=0;k<dim3;k++){
				array[i][j][k] = (int *****)malloc(sizeof(int****) * dim4);
				check_mem(array[i][j][k]);
				for(m=0;m<dim4;m++){
					array[i][j][k][m] = (int ****)malloc(sizeof(int***) * dim5);
					check_mem(array[i][j][k][m]);
					for(n=0;n<dim5;n++){					
						array[i][j][k][m][n] = (int ***)malloc(sizeof(int**) * dim6);
						check_mem(array[i][j][k][m][n]);
						for(p=0;p<dim6;p++){
							array[i][j][k][m][n][p] = (int **)malloc(sizeof(int*) * dim7);
							check_mem(array[i][j][k][m][n][p]);
							for(r=0;r<dim7;r++){
								array[i][j][k][m][n][p][r] = (int *)malloc(sizeof(int) * dim8);
								check_mem(array[i][j][k][m][n][p][r]);
								for(s=0;s<dim8;s++){
									array[i][j][k][m][n][p][r][s] = 0;
								}
							}
						}
					}
				}
			}
		}
	}
	return array;
error:
	return NULL;
}

void covs_access_free_cov_array_given_dimensions(int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8, int ********arr){
	int i,j,k,m,n,p,r;
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			for(k=0;k<dim3;k++){
				for(m=0;m<dim4;m++){
					for(n=0;n<dim5;n++){					
						for(p=0;p<dim6;p++){
							for(r=0;r<dim7;r++){
								free(arr[i][j][k][m][n][p][r]);
							}
							free(arr[i][j][k][m][n][p]);
						}
						free(arr[i][j][k][m][n]);
					}
					free(arr[i][j][k][m]);
				}
				free(arr[i][j][k]);
			}
			free(arr[i][j]);
		}
		free(arr[i]);
	}
	free(arr);
	return;
}

void covs_access_free_prob_array_given_dimensions(int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8, long double ********arr){
	int i,j,k,m,n,p,r;
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			for(k=0;k<dim3;k++){
				for(m=0;m<dim4;m++){
					for(n=0;n<dim5;n++){					
						for(p=0;p<dim6;p++){
							for(r=0;r<dim7;r++){
								free(arr[i][j][k][m][n][p][r]);
							}
							free(arr[i][j][k][m][n][p]);
						}
						free(arr[i][j][k][m][n]);
					}
					free(arr[i][j][k][m]);
				}
				free(arr[i][j][k]);
			}
			free(arr[i][j]);
		}
		free(arr[i]);
	}
	free(arr);
	return;
}

int covs_access_write_covs_to_file(char *file_loc,int ********arr,int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8){
	assert(file_loc != NULL);
	assert(arr != NULL);
	int true_arr[dim1][dim2][dim3][dim4][dim5][dim6][dim7][dim8];
	int i,j,k,m,n,p,r,s;
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			for(k=0;k<dim3;k++){
				for(m=0;m<dim4;m++){
					for(n=0;n<dim5;n++){					
						for(p=0;p<dim6;p++){
							for(r=0;r<dim7;r++){
								for(s=0;s<dim8;s++){
									true_arr[i][j][k][m][n][p][r][s] = arr[i][j][k][m][n][p][r][s];
								}
							}
						}
					}
				}
			}
		}
	}
	FILE *file = fopen(file_loc,"wb");
	check(file,"Error opening file to write cov array: %s.",file_loc);
	int chk = fwrite(true_arr,sizeof(int),dim1*dim2*dim3*dim4*dim5*dim6*dim7*dim8,file);	
	check(chk==dim1*dim2*dim3*dim4*dim5*dim6*dim7*dim8,"Error writing cov array to file.");
	fflush(file);
	fclose(file);
	return 0;
error:
	if(file) fclose(file);
	return -1;
}

int ********covs_access_read_covs_from_file(char *file_loc,int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8){
	assert(file_loc != NULL);
	//Use a real array rather than pointers for reading...
	int true_arr[dim1][dim2][dim3][dim4][dim5][dim6][dim7][dim8];
	FILE *file = fopen(file_loc,"rb");
	check(file,"Error opening file to read cov array: %s.",file_loc);
	int chk = fread(true_arr,sizeof(int),dim1*dim2*dim3*dim4*dim5*dim6*dim7*dim8,file);	
	check(chk==dim1*dim2*dim3*dim4*dim5*dim6*dim7*dim8,"Error reading cov array from file.");
	fclose(file);
	//Actual return array
	int ********arr = covs_access_generate_cov_array_given_dimensions(dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8);
	//copy values to array pointers.
	int i,j,k,m,n,p,r,s;
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			for(k=0;k<dim3;k++){
				for(m=0;m<dim4;m++){
					for(n=0;n<dim5;n++){					
						for(p=0;p<dim6;p++){
							for(r=0;r<dim7;r++){
								for(s=0;s<dim8;s++){
									arr[i][j][k][m][n][p][r][s] = true_arr[i][j][k][m][n][p][r][s];
								}
							}
						}
					}
				}
			}
		}
	}
	return arr;
error:
	if(file) fclose(file);
	return NULL;
}

int cov_access_compare_two_cov_arrays(int ********first,int ********second,int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8){
	int i,j,k,m,n,p,r,s;
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			for(k=0;k<dim3;k++){
				for(m=0;m<dim4;m++){
					for(n=0;n<dim5;n++){					
						for(p=0;p<dim6;p++){
							for(r=0;r<dim7;r++){
								for(s=0;s<dim8;s++){
									if(first[i][j][k][m][n][p][r][s] != second[i][j][k][m][n][p][r][s]){
										return -1;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

int cov_access_compare_two_prob_arrays(long double ********first,long double ********second,int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8){
	int i,j,k,m,n,p,r,s;
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			for(k=0;k<dim3;k++){
				for(m=0;m<dim4;m++){
					for(n=0;n<dim5;n++){					
						for(p=0;p<dim6;p++){
							for(r=0;r<dim7;r++){
								for(s=0;s<dim8;s++){
									if(abs(first[i][j][k][m][n][p][r][s] - second[i][j][k][m][n][p][r][s])>0.00001){
										return -1;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return 0;
}


void cov_access_print_cov_array(int ********arr,int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8){
	int i,j,k,m,n,p,r,s;
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			for(k=0;k<dim3;k++){
				for(m=0;m<dim4;m++){
					for(n=0;n<dim5;n++){					
						for(p=0;p<dim6;p++){
							for(r=0;r<dim7;r++){
								for(s=0;s<dim8;s++){
									printf("%d,%d,%d,%d,%d,%d,%d,%d - %d\n",i,j,k,m,n,p,r,s,arr[i][j][k][m][n][p][r][s]);
								}
							}
						}
					}
				}
			}
		}
	}
	printf("\n");
	return;
}	

void covs_access_merge_count_arrays(int ********arr_1, int ********arr_2, int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8){
	int i,j,k,m,n,p,r,s;
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			for(k=0;k<dim3;k++){
				for(m=0;m<dim4;m++){
					for(n=0;n<dim5;n++){					
						for(p=0;p<dim6;p++){
							for(r=0;r<dim7;r++){
								for(s=0;s<dim8;s++){
									arr_1[i][j][k][m][n][p][r][s] += arr_2[i][j][k][m][n][p][r][s];
								}
							}
						}
					}
				}
			}
		}
	}
	return;
}

long double ********covs_access_generate_probability_array(int ********counts,int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size){
	long double ********array = (long double ********)malloc(sizeof(long double *******) * rd_order_size);
	check_mem(array);
	int i,j,k,m,n,p,r,s;
	for(i=0;i<rd_order_size;i++){
		array[i] = (long double *******)malloc(sizeof(long double******) * strand_size);
		check_mem(array[i]);
		for(j=0;j<strand_size;j++){
			array[i][j] = (long double ******)malloc(sizeof(long double*****) * lane_size);
			check_mem(array[i][j]);
			for(k=0;k<lane_size;k++){
				array[i][j][k] = (long double *****)malloc(sizeof(long double****) * rpos_size);
				check_mem(array[i][j][k]);
				for(m=0;m<rpos_size;m++){
					array[i][j][k][m] = (long double ****)malloc(sizeof(long double***) * mq_size);
					check_mem(array[i][j][k][m]);
					for(n=0;n<mq_size;n++){					
						array[i][j][k][m][n] = (long double ***)malloc(sizeof(long double**) * bq_size);
						check_mem(array[i][j][k][m][n]);
						for(p=0;p<bq_size;p++){
							array[i][j][k][m][n][p] = (long double **)malloc(sizeof(long double*) * ref_b_size);
							check_mem(array[i][j][k][m][n][p]);
							for(r=0;r<ref_b_size;r++){
								array[i][j][k][m][n][p][r] = (long double *)malloc(sizeof(long double) * called_b_size);
								check_mem(array[i][j][k][m][n][p][r]);
								long sum_all = 0;
								int zeroCount = 0;
								//Iterate through the called bases.								
								for(s=0;s<called_b_size;s++){
									if(counts[i][j][k][m][n][p][r][s] == 0){
										counts[i][j][k][m][n][p][r][s] = 1;
										zeroCount++;
									}
									sum_all += counts[i][j][k][m][n][p][r][s];
									//array[i][j][k][m][n][p][r][s] = 0;
								}
								
								if(zeroCount < 4){
									for(s=0;s<called_b_size;s++){
										long double result = logl((((long double) counts[i][j][k][m][n][p][r][s]) / ((long double) sum_all)));
										array[i][j][k][m][n][p][r][s] = result;
									}
								}else{
									for(s=0;s<called_b_size;s++){
										//double result = log((double) counts[i][j][k][m][n][p][r][s] / (double) sum_all);
										//array[i][j][k][m][n][p][r][s] = result;
										long sump = 0;
										//If we had 4 zeros we iterate through each  readorder, base quality and read position so we have a pseudo count
										int read_ord = 0;
										int cnt = 0;
										for(read_ord = 0; read_ord<rd_order_size; read_ord++ ){
											int read_pos = 0;
											for(read_pos=0; read_pos<rpos_size; read_pos++){	
												int base_q = 0;
												for(base_q = 0; base_q<bq_size; base_q++){
													cnt++;
													if(counts[read_ord][j][k][read_pos][n][base_q][r][s] == 0){
														counts[read_ord][j][k][read_pos][n][base_q][r][s] = 1;
													}
													sump += counts[read_ord][j][k][read_pos][n][base_q][r][s];
												}
											}
										}
										counts[i][j][k][m][n][p][r][s] = sump;
									}		
									long sum = 0;							
									for(s=0;s<called_b_size;s++){
										sum += counts[i][j][k][m][n][p][r][s];
									}
									for(s=0;s<called_b_size;s++){
										long double result = logl((((long double) counts[i][j][k][m][n][p][r][s]) / ((long double) sum)));
										array[i][j][k][m][n][p][r][s] = result;
									}
								}									
							}
						}
					}
				}
			}
		}
	}
	return array;
error:
	return NULL;
}

void cov_access_print_prob_array(long double ********arr,int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8){
	int i,j,k,m,n,p,r,s;
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			for(k=0;k<dim3;k++){
				for(m=0;m<dim4;m++){
					for(n=0;n<dim5;n++){					
						for(p=0;p<dim6;p++){
							for(r=0;r<dim7;r++){
								for(s=0;s<dim8;s++){
									printf("%d,%d,%d,%d,%d,%d,%d,%d,%10.5Le\n",i,j,k,m,n,p,r,s,arr[i][j][k][m][n][p][r][s]);
								}
							}
						}
					}
				}
			}
		}
	}
	printf("\n");
	return;
}

void cov_access_print_cov_and_prob_array(int ********arr,long double ********arr_2,int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8){
	int i,j,k,m,n,p,r,s;
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			for(k=0;k<dim3;k++){
				for(m=0;m<dim4;m++){
					for(n=0;n<dim5;n++){					
						for(p=0;p<dim6;p++){
							for(r=0;r<dim7;r++){
								for(s=0;s<dim8;s++){
									printf("%d,%d,%d,%d,%d,%d,%d,%d\t%d\t%10.5Le\n",i,j,k,m,n,p,r,s,arr[i][j][k][m][n][p][r][s],arr_2[i][j][k][m][n][p][r][s]);
								}
							}
						}
					}
				}
			}
		}
	}
	printf("\n");
	return;
}

int covs_access_write_probs_to_file(char *file_loc,long double ********arr,int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8){
	assert(file_loc != NULL);
	assert(arr != NULL);
	long double true_arr[dim1][dim2][dim3][dim4][dim5][dim6][dim7][dim8];
	int i,j,k,m,n,p,r,s;
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			for(k=0;k<dim3;k++){
				for(m=0;m<dim4;m++){
					for(n=0;n<dim5;n++){					
						for(p=0;p<dim6;p++){
							for(r=0;r<dim7;r++){
								for(s=0;s<dim8;s++){
									true_arr[i][j][k][m][n][p][r][s] = arr[i][j][k][m][n][p][r][s];
								}
							}
						}
					}
				}
			}
		}
	}
	FILE *file = fopen(file_loc,"wb");
	check(file,"Error opening file to write cov array: %s.",file_loc);
	int chk = fwrite(true_arr,sizeof(long double),dim1*dim2*dim3*dim4*dim5*dim6*dim7*dim8,file);	
	check(chk==dim1*dim2*dim3*dim4*dim5*dim6*dim7*dim8,"Error writing prob array to file.");
	fflush(file);
	fclose(file);
	return 0;
error:
	if(file) fclose(file);
	return -1;
}

long double ********covs_access_generate_prob_array_given_dimensions(int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8){
	long double ********array = (long double ********)malloc(sizeof(long double *******) * dim1);
	check_mem(array);
	int i,j,k,m,n,p,r,s;
	for(i=0;i<dim1;i++){
		array[i] = (long double *******)malloc(sizeof(long double******) * dim2);
		check_mem(array[i]);
		for(j=0;j<dim2;j++){
			array[i][j] = (long double ******)malloc(sizeof(long double*****) * dim3);
			check_mem(array[i][j]);
			for(k=0;k<dim3;k++){
				array[i][j][k] = (long double *****)malloc(sizeof(long double****) * dim4);
				check_mem(array[i][j][k]);
				for(m=0;m<dim4;m++){
					array[i][j][k][m] = (long double ****)malloc(sizeof(long double***) * dim5);
					check_mem(array[i][j][k][m]);
					for(n=0;n<dim5;n++){					
						array[i][j][k][m][n] = (long double ***)malloc(sizeof(long double**) * dim6);
						check_mem(array[i][j][k][m][n]);
						for(p=0;p<dim6;p++){
							array[i][j][k][m][n][p] = (long double **)malloc(sizeof(long double*) * dim7);
							check_mem(array[i][j][k][m][n][p]);
							for(r=0;r<dim7;r++){
								array[i][j][k][m][n][p][r] = (long double *)malloc(sizeof(long double) * dim8);
								check_mem(array[i][j][k][m][n][p][r]);
								for(s=0;s<dim8;s++){
									array[i][j][k][m][n][p][r][s] = 0;
								}
							}
						}
					}
				}
			}
		}
	}
	return array;
error:
	return NULL;
}

long double  ********covs_access_read_probs_from_file(char *file_loc,int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8){
	assert(file_loc != NULL);
	//Use a real array rather than pointers for reading...
	long double true_arr[dim1][dim2][dim3][dim4][dim5][dim6][dim7][dim8];
	FILE *file = fopen(file_loc,"rb");
	check(file,"Error opening file to read cov array: %s.",file_loc);
	int chk = fread(true_arr,sizeof(long double),dim1*dim2*dim3*dim4*dim5*dim6*dim7*dim8,file);	
	check(chk==dim1*dim2*dim3*dim4*dim5*dim6*dim7*dim8,"Error reading cov array from file.");
	fclose(file);
	//Actual return array
	long double ********arr = covs_access_generate_prob_array_given_dimensions(dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8);
	//copy values to array pointers.
	int i,j,k,m,n,p,r,s;
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			for(k=0;k<dim3;k++){
				for(m=0;m<dim4;m++){
					for(n=0;n<dim5;n++){					
						for(p=0;p<dim6;p++){
							for(r=0;r<dim7;r++){
								for(s=0;s<dim8;s++){
									arr[i][j][k][m][n][p][r][s] = true_arr[i][j][k][m][n][p][r][s];
								}
							}
						}
					}
				}
			}
		}
	}
	return arr;
error:
	if(file) fclose(file);
	return NULL;
}