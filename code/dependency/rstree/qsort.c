
#include "qsort.h"

/*
  function : void q_swap(int *sorted_index, double *value, int i, int j)
  purpose  : to swap the two elements in both sorted index and value
             specified by the index i and index j 
  parameter: int *sorted_index - the array with two elements to be swapped
             double *value     - the array with two elements to be swapped
	     int i             - the index of swapping
	     int j             - the index of swapping
  return   : none
*/
void q_swap(int *sorted_index, double *value, int i, int j)
{ 
  int temp_index;
  double temp_value; 

  temp_value = value[i];  
  value[i] = value[j];
  value[j] = temp_value;
  
  temp_index = sorted_index[i];
  sorted_index[i] = sorted_index[j];
  sorted_index[j] = temp_index;
                
  return;
   
}


/*
  function : int partition(int *sorted_index, double *value, int start_index, int end_index)
  purpose  : to partition the sorted index and the value into two parts, which is used
             in quicksort
  parameter: int *sorted_index - the array of sorted index
             double *value     - the array of value
	     int start_index   - the index of the beginning of the array to be partitioned
	     int end_index     - the index of the end of the array to be partitioned
  return   : int - the index pointing to the partition index
*/
int partition(int *sorted_index, double *value, int start_index, int end_index)
{
  double pivot_value;  
  int i, j;
        
  pivot_value = value[start_index];
        
  i = start_index - 1;
  j = end_index + 1;
                
  while (1) {
   
        do {    
                j = j - 1;
  
        } while (value[j] < pivot_value);
 
        do {
                i = i + 1;

        } while (value[i] > pivot_value);

        if (i < j)
                q_swap(sorted_index, value, i, j);
        else
                return(j);
        
  }     


  // **** Ray have added the return value -1.
  return -1;

}



/* I think it is sorted in decreasing order */
/*
  function : void quicksort(int *sorted_index, double *value, int start_index, int end_index)
  purpose  : to sort the array value starting from start_index to end_index with the output
             in the sorted_index
  parameter: int *sorted_index - the sorted index of the array value
             double *value     - the array value to be sorted
	     int start_index   - the start index of the array to be sorted
	     int end_index     - the end index of the array to be sorted
  return   : none
*/
void quicksort(int *sorted_index, double *value, int start_index, int end_index)
{
  int pivot;

  if (start_index < end_index) {

        pivot = partition(sorted_index, value, start_index, end_index);

        quicksort(sorted_index, value, start_index, pivot);

        quicksort(sorted_index, value, pivot+1, end_index);

  }
 
  return;

}





/*
  function : void q_swap(int *sorted_index, double *value, int i, int j)
  purpose  : to swap the two elements in both sorted index and value
             specified by the index i and index j 
  parameter: int *sorted_index - the array with two elements to be swapped
             double *value     - the array with two elements to be swapped
	     int i             - the index of swapping
	     int j             - the index of swapping
  return   : none
*/
void q_swap_i(int *sorted_index, int *value, int i, int j)
{ 
  int temp_index;
  int temp_value; 

  temp_value = value[i];  
  value[i] = value[j];
  value[j] = temp_value;
  
  temp_index = sorted_index[i];
  sorted_index[i] = sorted_index[j];
  sorted_index[j] = temp_index;
                
  return;
   
}


/*
  function : int partition(int *sorted_index, double *value, int start_index, int end_index)
  purpose  : to partition the sorted index and the value into two parts, which is used
             in quicksort
  parameter: int *sorted_index - the array of sorted index
             double *value     - the array of value
	     int start_index   - the index of the beginning of the array to be partitioned
	     int end_index     - the index of the end of the array to be partitioned
  return   : int - the index pointing to the partition index
*/
int partition_i(int *sorted_index, int *value, int start_index, int end_index)
{
  int pivot_value;  
  int i, j;
        
  pivot_value = value[start_index];
        
  i = start_index - 1;
  j = end_index + 1;
                
  while (1) {
   
        do {    
                j = j - 1;
  
        } while (value[j] < pivot_value);
 
        do {
                i = i + 1;

        } while (value[i] > pivot_value);

        if (i < j)
                q_swap_i(sorted_index, value, i, j);
        else
                return(j);
        
  }     


  // **** Ray have added the return value -1.
  return -1;

}



/* I think it is sorted in decreasing order */
/*
  function : void quicksort(int *sorted_index, double *value, int start_index, int end_index)
  purpose  : to sort the array value starting from start_index to end_index with the output
             in the sorted_index
  parameter: int *sorted_index - the sorted index of the array value
             double *value     - the array value to be sorted
	     int start_index   - the start index of the array to be sorted
	     int end_index     - the end index of the array to be sorted
  return   : none
*/
void quicksort_i(int *sorted_index, int *value, int start_index, int end_index)
{
  int pivot;

  if (start_index < end_index) {

        pivot = partition_i(sorted_index, value, start_index, end_index);

        quicksort_i(sorted_index, value, start_index, pivot);

        quicksort_i(sorted_index, value, pivot+1, end_index);

  }
 
  return;

}








/*
  function : void q_swap(int *sorted_index, double *value, int i, int j)
  purpose  : to swap the two elements in both sorted index and value
             specified by the index i and index j 
  parameter: int *sorted_index - the array with two elements to be swapped
             double *value     - the array with two elements to be swapped
	     int i             - the index of swapping
	     int j             - the index of swapping
  return   : none
*/
void q_swap_i_indShort(short *sorted_index, int *value, int i, int j)
{ 
  short temp_index;
  int temp_value; 

  temp_value = value[i];  
  value[i] = value[j];
  value[j] = temp_value;
  
  temp_index = sorted_index[i];
  sorted_index[i] = sorted_index[j];
  sorted_index[j] = temp_index;
                
  return;
   
}


/*
  function : int partition(int *sorted_index, double *value, int start_index, int end_index)
  purpose  : to partition the sorted index and the value into two parts, which is used
             in quicksort
  parameter: int *sorted_index - the array of sorted index
             double *value     - the array of value
	     int start_index   - the index of the beginning of the array to be partitioned
	     int end_index     - the index of the end of the array to be partitioned
  return   : int - the index pointing to the partition index
*/
int partition_i_indShort(short *sorted_index, int *value, int start_index, int end_index)
{
  int pivot_value;  
  int i, j;
        
  pivot_value = value[start_index];
        
  i = start_index - 1;
  j = end_index + 1;
                
  while (1) {
   
        do {    
                j = j - 1;
  
        } while (value[j] < pivot_value);
 
        do {
                i = i + 1;

        } while (value[i] > pivot_value);

        if (i < j)
                q_swap_i_indShort(sorted_index, value, i, j);
        else
                return(j);
        
  }     


  // **** Ray have added the return value -1.
  return -1;

}



/* I think it is sorted in decreasing order */
/*
  function : void quicksort(int *sorted_index, double *value, int start_index, int end_index)
  purpose  : to sort the array value starting from start_index to end_index with the output
             in the sorted_index
  parameter: int *sorted_index - the sorted index of the array value
             double *value     - the array value to be sorted
	     int start_index   - the start index of the array to be sorted
	     int end_index     - the end index of the array to be sorted
  return   : none
*/
void quicksort_i_indShort(short *sorted_index, int *value, int start_index, int end_index)
{
  int pivot;

  if (start_index < end_index) {

        pivot = partition_i_indShort(sorted_index, value, start_index, end_index);

        quicksort_i_indShort(sorted_index, value, start_index, pivot);

        quicksort_i_indShort(sorted_index, value, pivot+1, end_index);

  }
 
  return;

}

/*

Value: 0 1 2 3 4 5
Index: 0 1 2 3 4 5
Sorting...
Value: 5 4 3 2 1 0
Index: 5 4 3 2 1 0


Value: 2 1 0 3 4 5
Index: 0 1 2 3 4 5
Sorting...
Value: 5 4 3 2 1 0
Index: 5 4 3 0 1 2

*/
void testQsort()
{
	int *value;
	int *index;

	int i;

	int noOfElt;

	noOfElt = 6;

	value = (int *) malloc(sizeof(int)*noOfElt);
	index = (int *) malloc(sizeof(int)*noOfElt);

	value[0] = 2;
	value[1] = 1;
	value[2] = 0;
	value[3] = 3;
	value[4] = 4;
	value[5] = 5;

	for (i = 0; i < noOfElt; i++)
	{
		index[i] = i;
	}

	printf("Value: ");
	for (i = 0; i < noOfElt; i++)
	{
		printf("%d ", value[i]);
	}
	printf("\n");

	printf("Index: ");
	for (i = 0; i < noOfElt; i++)
	{
		printf("%d ", index[i]);
	}
	printf("\n");

	printf("Sorting...\n");
	quicksort_i(index, value, 0, noOfElt-1);

	printf("Value: ");
	for (i = 0; i < noOfElt; i++)
	{
		printf("%d ", value[i]);
	}
	printf("\n");

	printf("Index: ");
	for (i = 0; i < noOfElt; i++)
	{
		printf("%d ", index[i]);
	}
	printf("\n");

	//{ int temp = 666023680;//2147483647;
	//printf("~~~~~~~~~~~ %d ~~~~~~~~~~~~\n", temp);
	//}

	free(value);
	free(index);
}
