

#ifndef QSORT_H
#define QSORT_H

#include <stdio.h>
#include <stdlib.h>

// decreasing order
void q_swap(int *sorted_index, double *value, int i, int j);
int partition(int *sorted_index, double *value, int start_index, int end_index);
void quicksort(int *sorted_index, double *value, int start_index, int end_index);

void q_swap_i(int *sorted_index, int *value, int i, int j);
int partition_i(int *sorted_index, int *value, int start_index, int end_index);
void quicksort_i(int *sorted_index, int *value, int start_index, int end_index);

void q_swap_i_indShort(short *sorted_index, int *value, int i, int j);
int partition_i_indShort(short *sorted_index, int *value, int start_index, int end_index);
void quicksort_i_indShort(short *sorted_index, int *value, int start_index, int end_index);

void testQsort();


#endif
