#ifndef __RSTREE_H_
#define __RSTREE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "qsort.h"
#include "donut.h"

#define DEBUG 1


#define CONFIG_FILE	"rtree.config"
#define SAVE_RTREE_FILE "save_rtree_file"

#define FALSE    	0
#define TRUE     	1

#define RANGE_SEARCH 0
#define kNN_SEARCH 1
#define CHOICE kNN_SEARCH

#define ASC_NUM  	48
#define NO_ID	 	-1
#define FOUND		1
#define NOT_FOUND 	0

#define ROOT  0   
#define LEAF  1
#define NODE  2

//#define INFINITY  FLT_MAX
#ifndef INFINITY
#define INFINITY  1E+37 //FLT_MAX
#endif
#define UNDEFINED -3  // for id of entries in PR

// #define R_FLOAT
#define R_TYPE int
#define R_LONG_TYPE long long
// #define R_TYPE float
// #define R_LONG_TYPE double

#define FLOAT_ZERO 0.001
#define INT_ZERO 100

/* Global variable ******************
m: min. number entries of each node;
M: max. number entries of each node;
dim: dimension of the rtree data.
*************************************/

typedef struct rtree_info_s
{
	int m, M, dim, reinsert_p, no_histogram;
	int extra_level;
} rtree_info;

typedef struct node { 
	R_TYPE *a;
	R_TYPE *b;
	int id;
	int attribute;
	int vacancy;
	struct node *parent;
	struct node **ptr; 
}   node_type;
	
typedef struct NN {	
	R_LONG_TYPE dist;
	int oid;
	struct node *pointer; 
	int level;
	struct NN *next; 
} NN_type;
		
typedef struct BranchArray {
	R_LONG_TYPE min;
	R_LONG_TYPE max;
	node_type *node;
} ABL;

typedef struct RangeReturn {
	R_LONG_TYPE dist;
	int oid;
	struct node *pointer;
	struct RangeReturn *prev;
	struct RangeReturn *next;
} RangeReturn_type;

typedef struct BundleReturn {
	int oid;
	int qid;
	struct node *pointer;
	struct BundleReturn *prev;
	struct BundleReturn *next;
} BundleReturn_type;

typedef struct NodeReturn {
	struct node *pointer;
	struct NodeReturn *prev;
	struct NodeReturn *next;
} NodeReturn_type;


typedef struct config { 
	int dim;
	int m;
	int M;
	int reinsert_p;
	int no_histogram;
	//char nodefile[FILENAME_MAX];
	//char rootfile[FILENAME_MAX];
	char queryfile[FILENAME_MAX];
	char positionfile[FILENAME_MAX]; 
	char save_tree_file[FILENAME_MAX]; 
}   config_type;



struct setNode_s
{
	int noOfNode;
	node_type **elt;
	int *level;
};

typedef struct setNode_s setNode;


void overflow(node_type *over_node, int over_level, int old_level, node_type
			  *extra_node, node_type *root, rtree_info *aInfo);

void initialize(config_type *config, char *configFile, rtree_info *aInfo);
void tree_node_allocate(node_type **node, rtree_info *aInfo);
void NN_freeChain(NN_type *aNN);

void build_tree(node_type **root, R_TYPE **data, int no_data, rtree_info *aInfo);
void k_NN_search(node_type *root, R_TYPE *query, int k, NN_type **returnResult, rtree_info *aInfo);

// rad, min, max are all squared
void k_NN_search_sphere(node_type *root, R_TYPE *query, int k, NN_type **returnResult, rtree_info *aInfo, R_LONG_TYPE rad);
void sphere_search(node_type *root, R_TYPE *query, R_LONG_TYPE min, R_LONG_TYPE max, RangeReturn_type **returnResult, rtree_info *r_info);

void NN_range_search(node_type *root, R_TYPE *query, RangeReturn_type **rr_ret, rtree_info *aInfo, R_LONG_TYPE rad, double span, int excl_list[], int excl_num);

void k_Donut_NN_search(node_type *root, R_TYPE *query, int k, NN_type **returnResult, rtree_info *aInfo);

// return total no. nodes accessed:
int rectangle_search(node_type *root, int no_query, R_TYPE **query_min, R_TYPE **query_max, RangeReturn_type **returnResult, rtree_info *aInfo);
// return total no. nodes accessed:
int rectangle_search_bundle(node_type *root, int no_query, R_TYPE **query, float error, float *min, float *max, BundleReturn_type **returnResult, rtree_info *aInfo);
// not used:
// int rectangle_search_bundle(node_type *root, float *min, float *max, node_type **returnResult, rtree_info *aInfo);

void save_rtree(node_type *root, const char save_tree_file[], rtree_info *aInfo);
void save_rtree_comb(node_type *root, FILE *fp, rtree_info *aInfo);
void read_rtree(node_type **root, const char save_tree_file[], rtree_info *aInfo);
void read_rtree_fp(node_type **root, FILE *fp, rtree_info *aInfo);
void free_tree(node_type *node, rtree_info *aInfo);

#endif