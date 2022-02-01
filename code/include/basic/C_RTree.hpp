#ifndef __C_RTREE_H_
#define __C_RTREE_H_

#include <string>
#include <fstream>
#include <vector>
#include <cstdio>
#include <unordered_set>
#include <stack>

extern "C" {
    #include "rtree.h"
}

#include "TriMesh.h"

class C_RTree
{
private:
	static const int RSTREE_SCALE = 1e5;

	static bool s_info_created;
	static rtree_info s_info;
	static std::string s_info_filename;

	static void read_rstree_info();

	node_type* m_root;

	R_TYPE* convert_pt(const trimesh::point& p) const;
	
	R_TYPE* convert_donut(const trimesh::point& o, float r, const trimesh::point& n) const;

public:
	C_RTree() {
		// if (!s_info_created) {
		// 	read_rstree_info();
		// 	s_info_created = true;
		// }
	}
	virtual ~C_RTree() {
		if (m_root) {
			// std::cout << "Free tree at " << m_root << std::endl;
		    free_tree(m_root, &s_info);
		}
	}

	void build(const trimesh::TriMesh* mesh);

	void read(FILE *fp);

	float knn_donut(const trimesh::point& o, float r, const trimesh::point& n, int k, int* ret, const std::unordered_set<int>& excl_id_set = {}) const;
	int    nn_donut(const trimesh::point& o, float r, const trimesh::point& n,                  const std::unordered_set<int>& excl_id_set = {}) const;

	float knn_sphere(const trimesh::point& p, float r, int k, int* ret, const std::unordered_set<int>& excl_id_set = {}) const;
	int    nn_sphere(const trimesh::point& p, float r,                  const std::unordered_set<int>& excl_id_set = {}) const;

	void range_sphere_min_max(const trimesh::point& p, float r_min, float r_max, std::vector<int>& ret, const std::unordered_set<int>& excl_id_set = {}) const;
	void range_sphere_dist_err(const trimesh::point& p, float dist, float err, std::vector<int>& ret, const std::unordered_set<int>& excl_id_set = {}) const;

	void nn_sphere_range(const trimesh::point& p, float sq_dist, float err, std::vector<int>& ret, int excl_id_list[] = {}, int excl_id_num = 0) const;

};


bool C_RTree::s_info_created = false;
rtree_info C_RTree::s_info = { 5, 10, 3, 7 }; // TODO: remove hard-coding
std::string C_RTree::s_info_filename = "rstree.pcd.config";

void C_RTree::read_rstree_info()
{
    std::ifstream ifs(s_info_filename);
    ifs >> s_info.m >> s_info.M >> s_info.dim >> s_info.reinsert_p;
    ifs.close();
}

R_TYPE* C_RTree::convert_pt(const trimesh::point& p) const
{
    R_TYPE* ret = (R_TYPE *) malloc(sizeof(R_TYPE) * 3);
    for (int i = 0; i < 3; i++) {
    	ret[i] = int(p[i] * RSTREE_SCALE);
    }
    return ret;
}

R_TYPE* C_RTree::convert_donut(const trimesh::point& o, float r, const trimesh::point& n) const
{
    R_TYPE* ret = (R_TYPE *) malloc(sizeof(R_TYPE) * 7);
    ret[0] = int(o[0] * RSTREE_SCALE);
    ret[1] = int(o[1] * RSTREE_SCALE);
    ret[2] = int(o[2] * RSTREE_SCALE);
    ret[3] = int(r * RSTREE_SCALE);
    ret[4] = int(n[0] * RSTREE_SCALE);
    ret[5] = int(n[1] * RSTREE_SCALE);
    ret[6] = int(n[2] * RSTREE_SCALE);
    return ret;
}


// public
void C_RTree::read(FILE *fp)
{
    read_rtree_fp(&m_root, fp, &s_info);
}

void C_RTree::build(const trimesh::TriMesh* mesh)
{
	int n = mesh->vertices.size();
	// R_TYPE** data = (R_TYPE **) malloc(sizeof(R_TYPE *) * n);
 //    for (int i = 0; i < n; i++) {
 //        data[i] = (R_TYPE *) malloc(sizeof(R_TYPE) * 3);
 //        for (int j = 0; j < 3; j++) {
 //            data[i][j] = (int) (mesh->vertices[i][j] * RSTREE_SCALE);
 //        }
 //    }

	R_TYPE** data = new R_TYPE *[n];
	for (int i = 0; i < n; i++) {
		data[i] = new R_TYPE[3];
		for (int j = 0; j < 3; j++) {
			data[i][j] = (int) (mesh->vertices[i][j] * RSTREE_SCALE);
		}
	}

	// build R-tree for mesh points
    // cout << "Start building R-tree for mesh #" << id << " of info #" << i << "..." << endl;
    build_tree(&m_root, data, n, &s_info);

    // // save R-tree
    // string realname_idx(realname + "." + to_string(i));
    // save_rtree(root, realname_idx.c_str(), &info[i]);
    // cout << "Save R-tree to " << realname_idx << endl;

	// for (int i = 0; i < n; i++) {
	// 	free(data[i]);
	// }
	// free(data);

	for (int i = 0; i < n; i++) {
		delete[] data[i];
	}
	delete[] data;
}

int C_RTree::nn_donut(const trimesh::point& o, float r, const trimesh::point& n, const std::unordered_set<int>& excl_id_set) const
{
	int ret;
	this->knn_donut(o, r, n, 1, &ret, excl_id_set);
	return ret;
}

float C_RTree::knn_donut(const trimesh::point& o, float r, const trimesh::point& n, int k, int* ret, const std::unordered_set<int>& excl_id_set) const
{
	auto query = convert_donut(o, r, n);
    NN_type *nn;
    int real_k = excl_id_set.size() + k;

    k_Donut_NN_search(m_root, query, real_k, &nn, &s_info);

    NN_type *nn_copy = nn;

    std::stack<NN_type*> s;
    for (int i = 0; i < real_k; i++) {
    	s.push(nn);
    	nn = nn->next;
    }

    float nn1_dist;
    int i = 0;
    while (i < k) {
    	auto top = s.top();
    	s.pop();
    	if (excl_id_set.find(top->oid) == excl_id_set.end()) { // not in the exclusive list
    		if (i == 0) {
    			nn1_dist = (float) top->dist / (float) RSTREE_SCALE;
    		}
    		ret[i] = top->oid;
    		i++;
    	}
    }

    NN_freeChain(nn_copy);

    return nn1_dist;
}

int C_RTree::nn_sphere(const trimesh::point& p, float r, const std::unordered_set<int>& excl_id_set) const
{
	int ret;
	this->knn_sphere(p, r, 1, &ret, excl_id_set);
	return ret;
}

float C_RTree::knn_sphere(const trimesh::point& p, float r, int k, int* ret, const std::unordered_set<int>& excl_id_set) const
{
	auto query = convert_pt(p);
    NN_type *nn;
    int real_k = excl_id_set.size() + k;

    if (r == 0.0) {
    	k_NN_search(m_root, query, real_k, &nn, &s_info);
    } else {
    	long long sq_r_long = (long long) (r * r * RSTREE_SCALE * RSTREE_SCALE);
    	k_NN_search_sphere(m_root, query, real_k, &nn, &s_info, sq_r_long);
    }

    NN_type *nn_copy = nn;

    std::stack<NN_type*> s;
    for (int i = 0; i < real_k; i++) {
    	s.push(nn);
    	nn = nn->next;
    }

    float nn1_dist;
    int i = 0;
    while (i < k) {
    	auto top = s.top();
    	s.pop();
    	if (excl_id_set.find(top->oid) == excl_id_set.end()) { // not in the exclusive list
    		if (i == 0) {
    			nn1_dist = (float) sqrt(float(top->dist) / float(RSTREE_SCALE) / float(RSTREE_SCALE));
    		}
    		ret[i] = top->oid;
    		i++;
    	}
    }

    NN_freeChain(nn_copy);

    return nn1_dist;
}

void C_RTree::range_sphere_min_max(const trimesh::point& p, float r_min, float r_max, std::vector<int>& ret, const std::unordered_set<int>& excl_id_set) const
{
	auto query = convert_pt(p);
	RangeReturn_type* rr = NULL;

	const long long sq_dist_min_long = (long long)(r_min * r_min * RSTREE_SCALE * RSTREE_SCALE);
	const long long sq_dist_max_long = (long long)(r_max * r_max * RSTREE_SCALE * RSTREE_SCALE);

	// cout << sq_dist_min_long << endl;
	// cout << sq_dist_max_long << endl;

	sphere_search(m_root, query, sq_dist_min_long, sq_dist_max_long, &rr, &s_info);

	while (rr) {
		if (excl_id_set.find(rr->oid) == excl_id_set.end()) { // not in the exclusive list
			// ret.push_back({ rr->oid, rr->dist });
			ret.push_back(rr->oid);
		}
		RangeReturn_type* prev = rr->prev;
		free(rr);
		rr = prev;
	}

	free(query);
}

void C_RTree::range_sphere_dist_err(const trimesh::point& p, float dist, float err, std::vector<int>& ret, const std::unordered_set<int>& excl_id_set) const
{
	range_sphere_min_max(p, std::max(0.0f, dist - err), dist + err, ret, excl_id_set);
}

void C_RTree::nn_sphere_range(const trimesh::point& p, float sq_dist, float err, std::vector<int>& ret, int excl_id_list[], int excl_id_num) const
{
	auto query = convert_pt(p);
	RangeReturn_type* rr = NULL;

	long long sq_dist_long;
	if (sq_dist == 0.0) {
		sq_dist_long = 0;
	} else {
		sq_dist_long = (long long) (sq_dist * RSTREE_SCALE * RSTREE_SCALE);
	}
	float span;
	if (err == 0.0) {
		span = 0.0;
	} else {
		span = err * 2 * RSTREE_SCALE;
	}

	NN_range_search(m_root, query, &rr, &s_info, sq_dist_long, span, excl_id_list, excl_id_num);

	while (rr) {
		// ret.push_back({ rr->oid, rr->dist });
		ret.push_back(rr->oid);
		RangeReturn_type* next = rr->next;
		free(rr);
		rr = next;
	}
}

#endif