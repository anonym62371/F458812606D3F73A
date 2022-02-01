
#include "rtree.h"

int E_dist_comp_count = 0;
int E_page_access_count = 0;
int E_represent_accuracy = 0;
int E_prune_ans_no = 0;
int E_overlap_yes = 0;
int E_no_answer = 0;



/*
void initialize(config_type *config, char *configFile)
{
FILE *fp;
int string_len;

 //fp = fopen(CONFIG_FILE, "r");
 fp = fopen(configFile, "r");
 
  
   fscanf(fp, "m=%d\n", &m);
   fscanf(fp, "M=%d\n", &M);
   fscanf(fp, "dim=%d\n", &dim);
   fscanf(fp, "reinsert_p=%d\n", &reinsert_p); 
   fscanf(fp, "no_histogram=%d\n",&no_histogram);
   
	
	 
	  fgets(config->positionfile, FILENAME_MAX, fp);
	  string_len = strlen(config->positionfile);
	  config->positionfile[string_len-1] = '\0';
	  
	   fgets(config->queryfile, FILENAME_MAX, fp);
	   string_len = strlen(config->queryfile);
	   config->queryfile[string_len-1] = '\0';
	   
		sprintf(config->save_tree_file, "%s.rstree", config->positionfile);
		}
		
		 
		  
		   void tree_node_allocate(node_type **node)
		   {
		   
			(*node) = (node_type *)malloc(sizeof(node_type));
			(*node)->a = (float *)malloc(sizeof(float) * dim);
			(*node)->b = (float *)malloc(sizeof(float) * dim);
			(*node)->ptr = (node_type **)malloc(sizeof(node_type *) * M);
			}
*/

int read_query(char *queryfile, R_TYPE ***query, int no_query, rtree_info *aInfo)
{
	FILE *fp_pos;  
	//int no_query;
	
	int i,j;
	
	//no_query=100;
	
	fp_pos = fopen(queryfile, "r");
	(*query) = (R_TYPE **)malloc(no_query * sizeof(R_TYPE*));
	for (i=0; i<no_query; i++)
		(*query)[i] = (R_TYPE *)malloc(aInfo->dim * sizeof(R_TYPE));
	
	for (i=0; i<no_query; i++)
		for (j=0; j<aInfo->dim; j++)
			fscanf(fp_pos, "%f", &((*query)[i][j]));
		fclose(fp_pos);
		
		return(no_query);
		
} /* read_query */

/* Distance Computation */
R_LONG_TYPE MINDIST(R_TYPE *P, R_TYPE *a, R_TYPE *b, rtree_info *aInfo)
{
	int i;

	R_LONG_TYPE diff;
	R_LONG_TYPE sum = 0.0;
	
	for(i=0 ;i<aInfo->dim ;i++)
	{
		if (P[i] > b[i])
		{
			//sum += pow(P[i] - b[i], 2.0);
			diff = P[i] - b[i];
			sum += diff*diff;
		}
		else if (P[i] < a[i])
		{
			//sum += pow(a[i] - P[i], 2.0);
			diff = a[i] - P[i];
			sum += diff*diff;
		}
	}
	
	return sum;     
}

R_LONG_TYPE MINDIST_Donut(struct Circle *cir, R_TYPE *a, R_TYPE *b, rtree_info *aInfo, double *time1, double *time2, double *time3)
{
	int i;

	struct Box box;
	for (i = 0; i < 3; i++) {
		box.a[i] = (double) a[i];
		box.b[i] = (double) b[i];
	}
	
	return (R_LONG_TYPE) min_dist_donut(&box, cir, time1, time2, time3);
}

R_LONG_TYPE MAXDIST(R_TYPE *P, R_TYPE *a, R_TYPE *b, rtree_info *aInfo)
{
	int i;
	R_TYPE mid;

	R_LONG_TYPE diff;
	R_LONG_TYPE sum = 0.0;

	for (i = 0; i < aInfo->dim; i++) {
#ifdef R_FLOAT
		mid = (a[i] + b[i]) * 0.5; // TODO: overflow?
#else
		mid = (a[i] + b[i] + 1) / 2;
#endif
		diff = (P[i] < mid) ? (b[i] - P[i]) : (P[i] - a[i]);
		sum += diff * diff;
	}

	return sum;
}

R_LONG_TYPE cal_Donut(node_type *node, struct Circle *cir, rtree_info *aInfo) {
	double x[3];
	int i;
	for (i = 0; i < 3; i++) {
		x[i] = (double) node->a[i];
	}
	return sqrt(donut_dist_sq(x, cir));
}

R_LONG_TYPE cal_Euclidean(node_type *node, R_TYPE *query, rtree_info *aInfo)
{
	int i;

	R_LONG_TYPE distance = 0.0;
	R_LONG_TYPE diff;
	
	for(i=0; i< aInfo->dim ;i++)
	{
		//distance += pow((node->a[i] - query[i]),(double)2.0);
		diff = node->a[i] - query[i];
		distance += diff*diff;
	}
	
	// return (sqrt(distance));
	return distance;
}

void RangeReturn_update(RangeReturn_type **RR, R_LONG_TYPE dist, node_type *node)
{
	RangeReturn_type* new_return = (RangeReturn_type *) malloc(sizeof(RangeReturn_type));
	new_return->dist = dist;
	new_return->oid = node->id;
	new_return->pointer = node;
	new_return->next = NULL;
	new_return->prev = *RR;

	if (*RR)
		(*RR)->next = new_return;
	
	*RR = new_return;
}


/***********************************/
/* rectangle_search_node():        */
/* search query points on the tree */
/***********************************/

void rectangle_search_node(node_type *curr_node, R_TYPE *query_min, R_TYPE *query_max, rtree_info *aInfo, RangeReturn_type **RR, int *page_accessed)
{
	int i, j, stop, flag;
	int query_dim;              
	
	query_dim = aInfo->dim;

	(*page_accessed)++;
	
	/* Search leaf node */
	if (curr_node->attribute == LEAF)
	{	
		RangeReturn_update(RR, 0, curr_node);	
	}
	
	stop = aInfo->M - curr_node->vacancy;
	for (i = 0; i < stop; i++)
	{
        flag = TRUE;
		
        /* search subtree */
        for (j = 0; j < query_dim; j++)
        {
			// if for any dimension, the current node box is not intersecting the query range,
			// prune this node
			if (curr_node->ptr[i]->a[j] > query_max[j] ||
				curr_node->ptr[i]->b[j] < query_min[j])
			{
				flag = FALSE;
				break;
			}
        }
		
        /* search the node which contains the query */
        if (flag == TRUE)
        {
			rectangle_search_node(curr_node->ptr[i], query_min, query_max, aInfo, RR, page_accessed);
        }
	}
}
/* rectangle_search_node */



/**********************************/
/* rectangle_search():            */
/* prepare to search query points */
/**********************************/

int rectangle_search(node_type *root, int no_query, R_TYPE **query_min, R_TYPE **query_max, RangeReturn_type **returnResult, rtree_info *aInfo)
{
	int query_index;
	int page_accessed = 0;
	
	/* start search data points by invoking rectangle_search_node() */
	for (query_index = 0; query_index < no_query; query_index++)
	{
        rectangle_search_node(root, query_min[query_index], query_max[query_index], aInfo, returnResult, &page_accessed);
	}

	return page_accessed;
} /* rectangle_search */


void BundleReturn_update(BundleReturn_type **BR, node_type *node, int qid)
{
	BundleReturn_type* new_return = (BundleReturn_type *) malloc(sizeof(BundleReturn_type));
	new_return->oid = node->id;
	new_return->qid = qid;
	new_return->pointer = node;
	new_return->next = NULL;
	new_return->prev = *BR;

	if (*BR)
		(*BR)->next = new_return;
	
	*BR = new_return;
}


void rectangle_search_node_bundle(node_type *curr_node, int no_query, R_TYPE **query, float error, float *min, float *max, rtree_info *aInfo, BundleReturn_type **BR, int *page_accessed)
{
	int i, j, k, stop, flag, flag_pr;
	R_TYPE *ind_query;
	int query_dim;            
	
	query_dim = aInfo->dim;

	(*page_accessed)++;
	
	/* Search leaf node */
	if (curr_node->attribute == LEAF)
	{
		for (i = 0; i < no_query; i++)
		{
			ind_query = query[i];
			flag = TRUE;
			for (j = 0; j < query_dim; j++)
			{
				if (curr_node->a[j] > (ind_query[j] + error) ||
					curr_node->a[j] < (ind_query[j] - error))
				{
					flag = FALSE;
					break;
				}
			}

			if (flag == TRUE)
			{
				BundleReturn_update(BR, curr_node, i);
			}
		}
	}
	
	stop = aInfo->M - curr_node->vacancy;
	for (i = 0; i < stop; i++)
	{
        flag = TRUE;
		
        /* search subtree */
        for (j = 0; j < query_dim; j++)
        {
			// if for any dimension, the current node box is not intersecting the ``big'' query range,
			// prune this node directly
			if (curr_node->ptr[i]->a[j] > max[j] ||
				curr_node->ptr[i]->b[j] < min[j])
			{
				flag = FALSE;
				break;
			}
        }

        if (flag == FALSE)
        	continue;

        flag = FALSE;

		for (j = 0; j < no_query; j++)
		{
			ind_query = query[j];
			flag_pr = TRUE;
			for (k = 0; k < query_dim; k++)
			{
				if (!(curr_node->a[k] > (ind_query[k] + error) ||
					  curr_node->b[k] < (ind_query[k] - error)))
				{
					flag_pr = FALSE;
					break;
				}
			}

			if (flag_pr == FALSE)
			{
				flag = TRUE;
				break;
			}
		}
		
        /* search the node which contains the query */
        if (flag == TRUE)
        {
        	rectangle_search_node_bundle(curr_node->ptr[i], no_query, query, error, min, max, aInfo, BR, page_accessed);
        }
	}
	
}


int rectangle_search_bundle(node_type *root, int no_query, R_TYPE **query, float error, float *min, float *max, BundleReturn_type **returnResult, rtree_info *aInfo)
{
	int page_accessed = 0;
	rectangle_search_node_bundle(root, no_query, query, error, min, max, aInfo, returnResult, &page_accessed);

	return page_accessed;
}


void read_inter_node(node_type *node, FILE *fp, rtree_info *aInfo)
{
	int i, count;
	
#ifdef R_FLOAT
	for (i = 0; i<aInfo->dim; i++)
        fscanf(fp, "%f\n", &((node->a)[i]));
	for (i = 0; i<aInfo->dim; i++)
        fscanf(fp, "%f\n", &((node->b)[i]));
#else
	for (i = 0; i<aInfo->dim; i++)
        fscanf(fp, "%d\n", &((node->a)[i]));
	for (i = 0; i<aInfo->dim; i++)
        fscanf(fp, "%d\n", &((node->b)[i]));
#endif
	
	fscanf(fp, "%d\n", &(node->attribute));
	if (node->attribute == LEAF) {
		fscanf(fp, "%d\n", &(node->id));
	}
	fscanf(fp, "%d\n", &(node->vacancy));
	
	/****************************************
	for (i = 0; i<aInfo->dim; i++)
	printf("%f\n", node->a[i]);
	
	 for (i = 0; i<aInfo->dim; i++)
	 printf("%f\n", node->b[i]);
	 
	  printf("%d\n", node->attribute);
	  if (node->attribute == LEAF) {
	  printf("%d\n", node->id);
	  }
	  printf("%d\n", node->vacancy);
	/****************************************/
	
	if (node->attribute != LEAF) {
		count = aInfo->M - node->vacancy;
		for (i=0; i<count; i++) {
			tree_node_allocate(&(node->ptr[i]), aInfo);
			read_inter_node(node->ptr[i], fp, aInfo);
		}
	}        
	
	return;
	
} 

void read_rtree_fp(node_type **root, FILE *fp, rtree_info *aInfo)
{
	// FILE *fp;
	
	tree_node_allocate(root, aInfo);
	
	//fp = fopen(SAVE_RTREE_FILE, "r");
	// fp = fopen(save_tree_file, "r");
	if (!fp) {
		printf("Can't open R-tree file\n");
		exit(EXIT_FAILURE);
	}

	read_inter_node(*root, fp, aInfo);
	
	// fclose(fp);
	
	return;

}


void read_rtree(node_type **root, const char save_tree_file[], rtree_info *aInfo)
{ 
	FILE *fp;
	
	
	tree_node_allocate(root, aInfo);
	
	//fp = fopen(SAVE_RTREE_FILE, "r");
	fp = fopen(save_tree_file, "r");
	if (!fp) {
		printf("Can't open file %s\n", save_tree_file);
		exit(EXIT_FAILURE);
	}

	read_inter_node(*root, fp, aInfo);
	
	fclose(fp);
	
	return;
	
}

void NN_freeChain(NN_type *NN)
{
	if (NN->next != NULL)
	{
		NN_freeChain(NN->next);
	}
	free(NN);
}

void NN_update(NN_type **NN, R_LONG_TYPE dist, node_type *node, int level, int k)
{
	int i=0;
	int l;
	int count;
	NN_type *nextNN;
	NN_type *oldNN;
	int noOfNode;
	
	if ((*NN)->dist == dist)
	{
		nextNN = (NN_type *) malloc(sizeof(NN_type));
		nextNN->oid = node->id;
		nextNN->dist = dist;
		nextNN->pointer = node;
		nextNN->level = level;
		nextNN->next = (*NN)->next;
		
		(*NN)->next = nextNN;
		
		return;
	}
	
	nextNN = (*NN);
	//while (i<k-1 && NN->next->dist >= dist)
	while ((nextNN->next != NULL) && (nextNN->next->dist > dist))
	{
		nextNN->dist = nextNN->next->dist;
		nextNN->level = nextNN->next->level;
		nextNN->oid = nextNN->next->oid;
		nextNN->pointer = nextNN->next->pointer;
		nextNN = nextNN->next;
		i++;
	}
	
	nextNN->dist = dist;
	nextNN->level = level;
	nextNN->oid = node->id;
	nextNN->pointer = node;
	
	
	// shrink
	noOfNode = 1;
	nextNN = (*NN);
	while (nextNN->next != NULL)
	{
		noOfNode++;
		nextNN = nextNN->next;
	}
	
	if (noOfNode != k)
	{
		//check
		l = noOfNode-k;
		
		// find l-th
		count = 1;
		nextNN = (*NN);
		while (count < l)
		{
			count++;
			nextNN = nextNN->next;
		}
		
		if (nextNN->next != NULL)
		{
			if (nextNN->dist == nextNN->next->dist)
			{
				// do nothing
			}
			else
			{
				oldNN = (*NN);
				(*NN) = nextNN->next;
				
				nextNN->next = NULL;
				
				
				// free
				NN_freeChain(oldNN);
			}
		}
		else
		{
			// impossible
		}
	}
	
	
	return;
}


//static int compare(ABL *i, ABL *j) 
//static int compare(const void *i, const void *j) // this form gets rid of warning
static int compare(const ABL *i, const ABL *j) 
{
	// const ABL* i_ABL = (ABL *) i;
	// const ABL* j_ABL = (ABL *) j;
	// if (i_ABL->min > j_ABL->min)
	// 	return (1);
	// if (i_ABL->min < j_ABL->min)
	// 	return (-1);
	if (i->min > j->min)
		return (1);
	if (i->min < j->min)
		return (-1);
	return (0);
}

void gen_ABL(node_type *node, ABL branch[], R_TYPE *query, int total, rtree_info *aInfo)
{
	int i;
	for (i=0;i<total;i++)
	{
		branch[i].node=node->ptr[i];
		branch[i].min=MINDIST(query, node->ptr[i]->a, node->ptr[i]->b, aInfo);
	}
	qsort(branch,total,sizeof(struct BranchArray),  compare);
	return;
}

void gen_ABL_minmax(node_type *node, ABL branch[], R_TYPE *query, int total, rtree_info *aInfo)
{
	int i;
	for (i=0;i<total;i++)
	{
		branch[i].node=node->ptr[i];
		branch[i].min=MINDIST(query, node->ptr[i]->a, node->ptr[i]->b, aInfo);
		branch[i].max=MAXDIST(query, node->ptr[i]->a, node->ptr[i]->b, aInfo);
	}
	qsort(branch,total,sizeof(struct BranchArray),  compare);
	return;
}

void gen_ABL_minmax_nosort(node_type *node, ABL branch[], R_TYPE *query, int total, rtree_info *aInfo)
{
	int i;
	for (i=0;i<total;i++)
	{
		branch[i].node=node->ptr[i];
		branch[i].min=MINDIST(query, node->ptr[i]->a, node->ptr[i]->b, aInfo);
		branch[i].max=MAXDIST(query, node->ptr[i]->a, node->ptr[i]->b, aInfo);
	}
	return;
}

void gen_Donut_ABL_minmax(node_type *node, ABL branch[], struct Circle *cir, int total, rtree_info *aInfo, double *time1, double *time2, double *time3)
{
	int i;
	for (i=0;i<total;i++)
	{
		branch[i].node=node->ptr[i];
		branch[i].min=MINDIST_Donut(cir, node->ptr[i]->a, node->ptr[i]->b, aInfo, time1, time2, time3);
		branch[i].max=MAXDIST(cir->o, node->ptr[i]->a, node->ptr[i]->b, aInfo);
	}
	qsort(branch,total,sizeof(struct BranchArray),  compare);
	return;
}

void k_NN_NodeSearch(node_type *curr_node, R_TYPE *query, NN_type **NN, int k, int level, rtree_info *aInfo)
{
	
	int i, total;
	R_LONG_TYPE dist;
	ABL *branch;
	
	E_page_access_count++;

	//*** new
	if (curr_node->vacancy == aInfo->M)
	{
		return;
	}
	
	/* Please refer NN Queries paper */
	if (curr_node->ptr[0]->attribute == LEAF) 
	{
		total = aInfo->M - curr_node->vacancy;
		for (i=0;i<total;i++)
		{
			dist = cal_Euclidean(curr_node->ptr[i], query, aInfo);
			if (dist <= (*NN)->dist)
				NN_update(NN, dist, curr_node->ptr[i], level+1, k);
		}
	}
	else
	{
		/* Please refer SIGMOD record Sep. 1998 Vol. 27 No. 3 P.18 */
		total=aInfo->M-curr_node->vacancy;
		branch=(struct BranchArray *)malloc(total*sizeof(struct BranchArray)*total);
		gen_ABL(curr_node, branch, query, total, aInfo);
		for (i=0;i<total;i++)
		{
			if (branch[i].min > (*NN)->dist)
				break;
			else
				k_NN_NodeSearch(branch[i].node, query, NN, k, level+1, aInfo);
		}
		free(branch);
	}
	return;
}


void k_NN_NodeSearch_sphere(node_type *curr_node, R_TYPE *query, NN_type **NN, int k, int level, rtree_info *aInfo, R_LONG_TYPE rad)
{
	
	int i, total;
	R_LONG_TYPE dist;
	ABL *branch;
	
	E_page_access_count++;

	//*** new
	if (curr_node->vacancy == aInfo->M)
	{
		return;
	}
	
	/* Please refer NN Queries paper */
	if (curr_node->ptr[0]->attribute == LEAF) 
	{
		total = aInfo->M - curr_node->vacancy;
		for (i=0;i<total;i++)
		{
			dist = cal_Euclidean(curr_node->ptr[i], query, aInfo);
			if (dist <= (*NN)->dist && dist > rad) // udpate only when dist > rad
				NN_update(NN, dist, curr_node->ptr[i], level+1, k);
		}
	}
	else
	{
		/* Please refer SIGMOD record Sep. 1998 Vol. 27 No. 3 P.18 */
		total=aInfo->M-curr_node->vacancy;
		branch=(struct BranchArray *)malloc(total*sizeof(struct BranchArray)*total);
		gen_ABL_minmax(curr_node, branch, query, total, aInfo);
		for (i=0;i<total;i++)
		{
			if (branch[i].min > (*NN)->dist)
				break;
			else if (branch[i].max <= rad)
				continue;
			else
				k_NN_NodeSearch_sphere(branch[i].node, query, NN, k, level+1, aInfo, rad);
		}
		free(branch);
	}
	return;

}

void k_Donut_NN_NodeSearch(node_type *curr_node, struct Circle *cir, NN_type **NN, int k, int level, rtree_info *aInfo, int *num_cal1, int *num_cal2, double *time1, double *time2, double *time3)
{
	
	int i, total;
	R_LONG_TYPE donut_dist, sq_eucl_dist;
	ABL *branch;
	
	E_page_access_count++;

	//*** new
	if (curr_node->vacancy == aInfo->M)
	{
		// printf("vacancy\n");
		return;
	}
	
	/* Please refer NN Queries paper */
	if (curr_node->ptr[0]->attribute == LEAF) 
	{
		total = aInfo->M - curr_node->vacancy;
		for (i=0;i<total;i++)
		{
			*(num_cal1) += 1;
			donut_dist = cal_Donut(curr_node->ptr[i], cir, aInfo);
			sq_eucl_dist = cal_Euclidean(curr_node->ptr[i], cir->o, aInfo);
			// printf("Donut dist leaf: %ld at %d\n", donut_dist, curr_node->ptr[i]->id);
			if (donut_dist < (*NN)->dist && sq_eucl_dist > cir->sq_r) {
				NN_update(NN, donut_dist, curr_node->ptr[i], level+1, k);
				// printf("Donut min update to %ld\n", (*NN)->dist);
			}
		}
	}
	else
	{
		/* Please refer SIGMOD record Sep. 1998 Vol. 27 No. 3 P.18 */
		total=aInfo->M-curr_node->vacancy;
		branch=(struct BranchArray *)malloc(total*sizeof(struct BranchArray)*total);
		gen_Donut_ABL_minmax(curr_node, branch, cir, total, aInfo, time1, time2, time3);
		*(num_cal2) += total;
		for (i=0;i<total;i++)
		{
			if (branch[i].min >= (*NN)->dist)
				break;
			else if (branch[i].max <= cir->sq_r)
				continue;
			else
				k_Donut_NN_NodeSearch(branch[i].node, cir, NN, k, level+1, aInfo, num_cal1, num_cal2, time1, time2, time3);
		}
		free(branch);
	}
	return;
}

void k_NN_trim(NN_type *head, int k)
{
	NN_type *aNN;
	NN_type *prevNN;
	NN_type *nextNN;
	NN_type *diffNN;
	double minDist;
	
	int i;
	aNN = head;
	
	for (i = 0; i < k-1; i++)
	{
		aNN = aNN->next;
	}
	
	minDist = aNN->dist;
	
	prevNN = aNN;
	nextNN = prevNN->next;
	diffNN = NULL;
	while (nextNN != NULL)
	{
		if (nextNN->dist != minDist)
		{
			diffNN = nextNN;
			
			prevNN->next = NULL;
			break;
		}
		
		prevNN = nextNN;
		nextNN = prevNN->next;
	}
	
	if (diffNN != NULL)
	{
		// free the remaining
		NN_freeChain(diffNN);
	}
}

void k_NN_search(node_type *root, R_TYPE *query, int k, NN_type **returnResult, rtree_info *aInfo)
{
	int i;
	
	NN_type *NN, *head;
	
	if ((NN = (NN_type *)malloc(sizeof(NN_type))) == NULL)
        fprintf(stderr, "malloc error at k-NN_search 1\n");
	NN->oid = UNDEFINED;
	NN->level = -1;
	NN->dist = INFINITY;
	NN->pointer = NULL;
	head=NN;
	
	for (i=0; i<k-1; i++) {
		
        if ((NN->next = (NN_type *)malloc(sizeof(NN_type))) == NULL)
			fprintf(stderr, "malloc error at k-NN_search 1\n");
		
        NN->next->oid = UNDEFINED; 
		NN->next->level = -1;
        NN->next->dist = INFINITY;
		NN->next->pointer = NULL;
		NN->next->next = NULL;
		
        NN = NN->next;
	}
	NN->next = NULL;
	
	k_NN_NodeSearch(root, query, &head, k, 0, aInfo);
	//k_NN_trim(head, k);
	
	NN=head;
/*
	printf("\n");
	
	while (NN != NULL) {
		
        printf("ID: %d    Dist: %f\n", NN->oid, NN->dist);
		
        NN = NN->next;
		
	}
	printf("\n");
*/
	
	(*returnResult) = head;
	
	return;
	
}

/**
 * Give me the squred radius
 */
void k_NN_search_sphere(node_type *root, R_TYPE *query, int k, NN_type **returnResult, rtree_info *aInfo, R_LONG_TYPE rad)
{
	int i;
	
	NN_type *NN, *head;
	
	if ((NN = (NN_type *)malloc(sizeof(NN_type))) == NULL)
        fprintf(stderr, "malloc error at k-NN_search 1\n");
	NN->oid = UNDEFINED;
	NN->level = -1;
	NN->dist = INFINITY;
	NN->pointer = NULL;
	head=NN;
	
	for (i=0; i<k-1; i++) {
		
        if ((NN->next = (NN_type *)malloc(sizeof(NN_type))) == NULL)
			fprintf(stderr, "malloc error at k-NN_search 1\n");
		
        NN->next->oid = UNDEFINED; 
		NN->next->level = -1;
        NN->next->dist = INFINITY;
		NN->next->pointer = NULL;
		NN->next->next = NULL;
		
        NN = NN->next;
	}
	NN->next = NULL;
	
	k_NN_NodeSearch_sphere(root, query, &head, k, 0, aInfo, rad);
	//k_NN_trim(head, k);
	
	NN=head;
	
	(*returnResult) = head;
	
	return;
}

void k_Donut_NN_search(node_type *root, R_TYPE *query, int k, NN_type **returnResult, rtree_info *aInfo)
{
	int i;
	
	NN_type *NN, *head;
	
	if ((NN = (NN_type *)malloc(sizeof(NN_type))) == NULL)
        fprintf(stderr, "malloc error at k-NN_search 1\n");
	NN->oid = UNDEFINED;
	NN->level = -1;
	NN->dist = INFINITY;
	NN->pointer = NULL;
	head=NN;
	
	for (i=0; i<k-1; i++) {
		
        if ((NN->next = (NN_type *)malloc(sizeof(NN_type))) == NULL)
			fprintf(stderr, "malloc error at k-NN_search 1\n");
		
        NN->next->oid = UNDEFINED; 
		NN->next->level = -1;
        NN->next->dist = INFINITY;
		NN->next->pointer = NULL;
		NN->next->next = NULL;
		
        NN = NN->next;
	}
	NN->next = NULL;

	struct Circle query_cir;
	query_cir.o[0] = (double) query[0];
	query_cir.o[1] = (double) query[1];
	query_cir.o[2] = (double) query[2];
	query_cir.r    = (double) query[3];
	query_cir.n[0] = (double) query[4];
	query_cir.n[1] = (double) query[5];
	query_cir.n[2] = (double) query[6];
	// for (i=0;i<7;i++) {
	// 	printf("%d\n", query[i]);
	// }

	int num1 = 0;
	int num2 = 0;
	double time1 = 0.0;
	double time2 = 0.0;
	double time3 = 0.0;

	pre_compute(&query_cir);	
	k_Donut_NN_NodeSearch(root, &query_cir, &head, k, 0, aInfo, &num1, &num2, &time1, &time2, &time3);
	// printf("%ld at %d\n", head->dist, head->oid);
	// printf("cal leaf: %d, cal node: %d, time: %f, %f, %f\n", num1, num2, time1, time2, time3);
	
	//k_NN_trim(head, k);
	
	NN=head;
/*
	printf("\n");
	
	while (NN != NULL) {
		
        printf("ID: %d    Dist: %f\n", NN->oid, NN->dist);
		
        NN = NN->next;
		
	}
	printf("\n");
*/
	
	(*returnResult) = head;
	
	return;
	
}

void sphere_NodeSearch(node_type *curr_node, R_TYPE *query, R_LONG_TYPE min, R_LONG_TYPE max, RangeReturn_type **RR, rtree_info *aInfo)
{
	
	int i, total;
	R_LONG_TYPE dist;
	ABL *branch;
	
	E_page_access_count++;

	//*** new
	if (curr_node->vacancy == aInfo->M)
	{
		return;
	}
	
	if (curr_node->ptr[0]->attribute == LEAF) 
	{
		total = aInfo->M - curr_node->vacancy;
		for (i=0;i<total;i++)
		{
			dist = cal_Euclidean(curr_node->ptr[i], query, aInfo);
#ifdef R_FLOAT
			if (min - FLOAT_ZERO <= dist && dist <= max + FLOAT_ZERO)
#else
			if (min - INT_ZERO <= dist && dist <= max + INT_ZERO)
#endif
				RangeReturn_update(RR, dist, curr_node->ptr[i]);
		}
	}
	else
	{
		/* Please refer SIGMOD record Sep. 1998 Vol. 27 No. 3 P.18 */
		total=aInfo->M-curr_node->vacancy;
		branch=(struct BranchArray *)malloc(total*sizeof(struct BranchArray)*total);
		gen_ABL_minmax_nosort(curr_node, branch, query, total, aInfo);
		for (i=0;i<total;i++)
		{
#ifdef R_FLOAT
			if (branch[i].min > max + FLOAT_ZERO || branch[i].max < min - FLOAT_ZERO)
#else
			if (branch[i].min > max + INT_ZERO || branch[i].max < min - INT_ZERO)
#endif
				continue;
			else
				sphere_NodeSearch(branch[i].node, query, min, max, RR, aInfo);
		}
		free(branch);
	}
	return;
}

/**
 * Gives me the squred distance
 */
void sphere_search(node_type *root, R_TYPE *query, R_LONG_TYPE min, R_LONG_TYPE max, RangeReturn_type **returnResult, rtree_info *r_info)
{
	*returnResult = NULL;
	sphere_NodeSearch(root, query, min, max, returnResult, r_info);
}


int ll_length(RangeReturn_type **rr_ret)
{
	int len = 0;
	RangeReturn_type *itr = *rr_ret;
	while (itr) {
		len++;
		itr = itr->next;
	}
	return len;
}


void NN_range_update(RangeReturn_type **rr_ret, R_LONG_TYPE dist, node_type *node, R_LONG_TYPE *nn_min, R_LONG_TYPE *nn_max, double span)
{

	// TODO: print min
	// printf("Update: %ld\n", dist);

	if (dist < *nn_min)
	{
		// update nn_min & nn_max
		*nn_min = dist;
		double sqrt_max = sqrt((double) dist) + span;
		*nn_max = (R_LONG_TYPE) (sqrt_max * sqrt_max) + INT_ZERO;
		// printf("Updated min: %ld\n", *nn_min);
		// printf("Updated max: %ld\n", *nn_max);

		// delete following rr_ret nodes which > updated *nn_max
		RangeReturn_type *itr = *rr_ret;
		while (itr)
		{
			if (itr->dist > *nn_max)
				break;
			itr = itr->next;
		}
		if (itr)
		{
			if (itr->prev)
			{
				itr->prev->next = NULL;
			}
			else
			{
				*rr_ret = NULL;
			}
			// while (itr != NULL)
			// TODO: delete the following, but first not do it, test memory and time
		}
	}

	RangeReturn_type* new_return = (RangeReturn_type *) malloc(sizeof(RangeReturn_type));
	new_return->dist = dist;
	new_return->oid = node->id;
	new_return->pointer = node;
	new_return->next = NULL;
	new_return->prev = NULL;

	// insert new_return into a proper position
	RangeReturn_type *itr = *rr_ret;
	RangeReturn_type *itr_prev = NULL;
	while (itr) {
		if (itr->dist >= dist)
		{
			break;
		}
		itr_prev = itr;
		itr = itr->next;
	}

	if (itr) // has next, link new to next
	{
		new_return->next = itr;
		itr->prev = new_return;
	}
	if (itr_prev) // has prev, link new to prev
	{
		itr_prev->next = new_return;
		new_return->prev = itr_prev;
	}
	else // has no prev, meaning new is the first, link *rr_ret to new
	{
		*rr_ret = new_return;
	}

	// printf("Length: %d\n", ll_length(rr_ret));
}

int contain(int excl_list[], int excl_num, int id)
{
	int i;
	for (i = 0; i < excl_num; i++)
	{
		if (excl_list[i] == id)
		{
			return 1;
		}
	}
	return 0;
}


void NN_range_NodeSearch(node_type *curr_node, R_TYPE *query, RangeReturn_type **rr_ret, rtree_info *aInfo, R_LONG_TYPE *nn_min, R_LONG_TYPE *nn_max,
	R_LONG_TYPE rad, double span, int excl_list[], int excl_num)
{
	
	int i, total;
	R_LONG_TYPE dist;
	ABL *branch;
	
	E_page_access_count++;

	//*** new
	if (curr_node->vacancy == aInfo->M)
	{
		return;
	}
	
	/* Please refer NN Queries paper */
	if (curr_node->ptr[0]->attribute == LEAF) 
	{
		total = aInfo->M - curr_node->vacancy;
		for (i=0;i<total;i++)
		{
			dist = cal_Euclidean(curr_node->ptr[i], query, aInfo);
			if (dist <= *nn_max && dist > rad && !contain(excl_list, excl_num, curr_node->ptr[i]->id)) // udpate only when dist > rad
				NN_range_update(rr_ret, dist, curr_node->ptr[i], nn_min, nn_max, span);
		}
	}
	else
	{
		/* Please refer SIGMOD record Sep. 1998 Vol. 27 No. 3 P.18 */
		total=aInfo->M-curr_node->vacancy;
		branch=(struct BranchArray *)malloc(total*sizeof(struct BranchArray)*total);
		gen_ABL_minmax(curr_node, branch, query, total, aInfo);
		for (i=0;i<total;i++)
		{
			if (branch[i].min > *nn_max)
				break;
			else if (branch[i].max <= rad)
				continue;
			else
				NN_range_NodeSearch(branch[i].node, query, rr_ret, aInfo, nn_min, nn_max, rad, span, excl_list, excl_num);
		}
		free(branch);
	}
	return;

}


/**
 * Give me the squred distance
 */
void NN_range_search(node_type *root, R_TYPE *query, RangeReturn_type **rr_ret, rtree_info *aInfo, R_LONG_TYPE rad, double span, int excl_list[], int excl_num)
{
	*rr_ret = NULL;

	R_LONG_TYPE *nn_min = (R_LONG_TYPE *) malloc(sizeof(R_LONG_TYPE));
	R_LONG_TYPE *nn_max = (R_LONG_TYPE *) malloc(sizeof(R_LONG_TYPE));
	*nn_min = INFINITY;
	*nn_max = INFINITY;

	NN_range_NodeSearch(root, query, rr_ret, aInfo, nn_min, nn_max, rad, span, excl_list, excl_num);

	free(nn_min);
	free(nn_max);

}


void free_tree(node_type *node, rtree_info *aInfo)
{
	int i;
	
	if (!node) return;
	
	free(node->a);
	free(node->b);
	
	for (i=0; i < aInfo->M-node->vacancy; i++)
        free_tree(node->ptr[i], aInfo);
	
	free(node);
	
	return;
}


// int main2(int argc, char *argv[])
// {
// 	rtree_info aInfo;
// 	node_type *root;
// 	config_type config;
// 	int no_query;
// 	R_TYPE **query;
// 	float error;
	
// 	int i, nn;
// 	NN_type *result;  
// 	//char E_result_filename[FILENAME_MAX], dummy[80];
	
// #ifndef WIN32
// 	// for experiment
// 	float  userTime, sysTime;
// 	struct rusage myTime1, myTime2;
// #endif
// 	FILE *E_result_fp;            
// 	double E_represent_accuracy_sum, E_prune_accuracy_sum;
// 	int E_page_access_count_sum;
	
// 	//////////////
// 	if (argc != 4) {
// 		printf("Usage: %s <config> <k> <#query> <result file>\n", argv[0]);
// 		exit(EXIT_FAILURE);
// 	}
	
// 	//strcpy(E_result_filename, argv[4]);
	
// 	//////////////
// 	initialize(&config, argv[1], &aInfo);
	
// 	printf("read_rtree\n");
// 	read_rtree(&root, config.save_tree_file, &aInfo);
	
// 	printf("read_query\n");
// 	no_query = read_query(config.queryfile, &query, atoi(argv[3]), &aInfo);
	
// 	printf("start searching\n");
	
// 	if (CHOICE==RANGE_SEARCH) {
		
// 		while (TRUE) {
			
// 			printf("Please input the error bound:");
// 			scanf("%f", &error);
// 			rectangle_search_tree(root, no_query, query, error, &aInfo);
			
// 		}
// 	}
	
// 	else if (CHOICE == kNN_SEARCH) {
		
// 		nn = atoi(argv[2]);
		
// 		E_page_access_count_sum = 0;
// 		E_represent_accuracy_sum = 0.0; 
// 		E_prune_accuracy_sum = 0.0;
// 		E_dist_comp_count = 0;
// #ifndef WIN32                         
// 		getrusage(RUSAGE_SELF,&myTime1);
// #endif
// 		for (i=0; i<no_query; i++) {
			
// 			printf("Query %d\n", i+1);
// 			E_page_access_count = 0;
// 			E_represent_accuracy = 0;
// 			E_prune_ans_no = 0;
// 			E_overlap_yes = 0;
			
// 			k_NN_search(root, query[i], nn, &result, &aInfo);
			
// 			E_page_access_count_sum += E_page_access_count;
// 			E_represent_accuracy_sum += ((double)E_represent_accuracy)/((double)E_page_access_count);
// 			E_prune_accuracy_sum += ((double)E_overlap_yes)/((double)E_prune_ans_no);
			
// 			printf("---------------------------\n");
// 		}
// #ifndef WIN32                                
// 		getrusage(RUSAGE_SELF,&myTime2);
		
//         userTime =
// 			((float) (myTime2.ru_utime.tv_sec  - myTime1.ru_utime.tv_sec)) +
// 			((float) (myTime2.ru_utime.tv_usec - myTime1.ru_utime.tv_usec)) * 1e-6;
//         sysTime =
// 			((float) (myTime2.ru_stime.tv_sec  - myTime1.ru_stime.tv_sec)) +
// 			((float) (myTime2.ru_stime.tv_usec - myTime1.ru_stime.tv_usec)) * 1e-6;
// #endif
		
//         //E_result_fp = fopen(E_result_filename, "w");
//         E_result_fp = stdout;
		
// 		fprintf(E_result_fp, "Number of query is %d\n", no_query);
// 		fprintf(E_result_fp, "average no. of page accesses = %f\n", (double)E_page_access_count_sum/(double)no_query);
		
// #ifndef WIN32
//         fprintf(E_result_fp, "User time : %f seconds\n",userTime);
//         fprintf(E_result_fp, "System time : %f seconds\n",sysTime);
//         fprintf(E_result_fp, "Total time : %f seconds\n",userTime+sysTime);
// #endif
		
// 		fprintf(E_result_fp, "\n");
// 		fclose(E_result_fp);
		
// 	}
	
// 	printf("free_tree\n");
// 	free_tree(root, &aInfo);
	
// 	return(0);
// }
