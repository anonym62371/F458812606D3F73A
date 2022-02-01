#ifndef __C2O_H_
#define __C2O_H_

#include <iostream>
#include <cassert>
#include <unordered_set>

#include "TriMesh.h"

#include "C_RTree.hpp"
#include "CPP_RTree.hpp"

// #define DENSE_THRESH 1e5 // such that only the 1M and 10M will be marked dense

#ifdef _IDX3
	#ifdef _CLR
		#define INDEX_DIM 6
	#else
		#define INDEX_DIM 3
	#endif
#else
	#define INDEX_DIM 6
#endif

// num of rd-filtering parameters, currently we have 6 sides
// in future, this will be designed to include user-customized parameters, e.g., color
#ifdef _CLR
	#define RD_NUM 7
#else
	#define RD_NUM 6
#endif

typedef RTree<int, int, INDEX_DIM, double, 36> IndexTree;
typedef RTree<int, int, INDEX_DIM, double, 6> AuxTree;
typedef RTree<int, int, 3, double, 6> CollTree;

enum C2OVariantEnum {
	C2OV_DONUT, C2OV_3NN, C2OV_3NN_SIM
};

trimesh::point getNNPtOnDonut(const trimesh::point& x, const trimesh::point& m, const trimesh::vec& n, float r_donut)
{
	auto len2_n = len2(n);
    auto mx = x - m;
    auto u = m + n * dot(n, mx) / len2_n; // let y be the return, u is the point on n, such that ux // my (same direction)
    auto ux = x - u;
    auto len_ux = len(ux);
    if (len_ux != 0.0f) {
        return (m + ux * (r_donut / len_ux));
    } else {
        // it means that x is accidentally on n, then u and x are the same point
        // then we pick arbitrary point on donut
        int idx = -1;
        float val = 0;
        for (int i = 0; i < 3; i++) {
            if (n[i] != 0.0f) {
                idx = i;
                val = (n[(i+1)%3] + n[(i+2)%3]) / n[i];
                break;
            }
        }
        assert(idx >= 0);
        trimesh::vec naive_parallel;
        naive_parallel[idx] = val;
        naive_parallel[(idx+1)%3] = -1;
        naive_parallel[(idx+2)%3] = -1;
        return (m + naive_parallel * (r_donut / len_ux));
    }
}

float donutDist2(const trimesh::point& x, const trimesh::point& m, const trimesh::vec& n, float r_donut)
{
    float sq_a = dist2(x, m);
    float sq_n = len2(n);
    float lambda = (dot(m, n) - dot(x, n)) / sq_n;
    float sq_h = sq(lambda) * sq_n;

    return (sq_a + sq(r_donut) - 2.0 * r_donut * sqrt(sq_a - sq_h));
}

// a quick solution, applied for sparse dataset
int quickDonutNNSparse(const trimesh::point& m, const trimesh::vec& n, float r_donut, const trimesh::TriMesh* mesh, const C_RTree* rtree)
{
    std::vector<int> donut_nn_cand;
    int order = 1;
    while (donut_nn_cand.size() == 0 && order <= 10) {
    	rtree->range_sphere_min_max(m, r_donut * float(order), r_donut * float(order + 1), donut_nn_cand); // automatically exclude p and a
    	order++;
    }

    // std::cout << "Order: " << order << std::endl;
    // std::cout << donut_nn_cand.size() << std::endl;
    if (donut_nn_cand.size() <= 0) {
        return -1;
    }

    int donut_nn_i = -1;
    float sq_d_min = std::numeric_limits<float>::max();

    for (auto &i: donut_nn_cand) {
        if (i < 0 || i >= mesh->vertices.size()) {
            continue;
        }
        const trimesh::point& i_pt = mesh->vertices[i];

        // auto pi = *i_pt.pt - *p->pt;
        // if (cos_theta(n, &pi) > 0.5) {
        //  continue;
        // }

        const float sq_d = donutDist2(i_pt, m, n, r_donut);

        if (sq_d < sq_d_min) {
            sq_d_min = sq_d;
            donut_nn_i = i;
        }
    }

    return donut_nn_i;
}

int quickDonutNNDense(const trimesh::point& m, const trimesh::vec& n, float r_donut, const trimesh::TriMesh* mesh, const C_RTree* rtree)
{
    int donut_nn_i;
    float sq_d_min = std::numeric_limits<float>::max();
    float d_min;

    const float intv = r_donut * 0.02;
    for (int t = 0; t < 5; t++) {
        std::vector<int> donut_nn_cand;
        const float left = float(t) * intv;
        const float rght = float(t + 1) * intv;
        rtree->range_sphere_min_max(m, r_donut + left, r_donut + rght, donut_nn_cand); // automatically exclude p and a
        // cout << "Itr size: " << donut_nn_cand.size() << endl;

        for (auto &i: donut_nn_cand) {
            if (i < 0 || i >= mesh->vertices.size()) {
                continue;
            }
        	const trimesh::point& i_pt = mesh->vertices[i];

            // auto pi = *i_pt.pt - *p->pt;
            // if (cos_theta(n, &pi) > 0.5) {
            //  continue;
            // }

            const float sq_d = donutDist2(i_pt, m, n, r_donut);

            if (sq_d < sq_d_min) {
                sq_d_min = sq_d;
                donut_nn_i = i;
            }
        }

        d_min = sqrt(sq_d_min);
        // cout << "Itr d_min: " << d_min << endl;
        if (d_min <= rght) {
            return donut_nn_i;
        }
    }

    return donut_nn_i;
}

int quickDonutNN(const trimesh::point& m, const trimesh::vec& n, float r_donut, const trimesh::TriMesh* mesh, const C_RTree* rtree, bool denseFlag)
{
	if (denseFlag) {
		return quickDonutNNDense(m, n, r_donut, mesh, rtree);
	} else {
		return quickDonutNNSparse(m, n, r_donut, mesh, rtree);
	}
}

// rotate point a anti-clockwise, in the view point reversed to vector n
trimesh::point getPtBOnDonut(const trimesh::point& a, const trimesh::point& m, const trimesh::vec& n)
{
    auto ma = a - m; // r_donut is supposed to equal ||ma||
    auto e = m + ma / 3.0;
    auto eb_dir = cross(n, ma); // such that the direction is ensured
    auto eb = eb_dir * 0.942809041 / len(n); // such that the length of eb is ensured, to verify it, mb = ma = r_donut

    return (e + eb);
}

bool calDonutEntry(int ptIDs[4], float r, const trimesh::TriMesh* mesh, const C_RTree* rtree, bool denseFlag, bool verbose)
{
    const trimesh::point& p1 = mesh->vertices[ptIDs[0]];

    ptIDs[1] = rtree->nn_sphere(p1, r);
    if (ptIDs[1] < 0) {
    	return false;
    }

    const trimesh::point& p2 = mesh->vertices[ptIDs[1]];

    auto m = middle_pt(p1, p2);
    auto n = p1 - p2;

    const float r_donut = dist(p1, m) * SQRT_3;

    ptIDs[2] = quickDonutNN(m, n, r_donut, mesh, rtree, denseFlag);
    if (ptIDs[2] < 0) {
    	return false;
    }

    const trimesh::point& p3 = mesh->vertices[ptIDs[2]];

    auto a = getNNPtOnDonut(p3, m, n, r_donut);
    auto b = getPtBOnDonut(a, m, n);

    ptIDs[3] = rtree->nn_sphere(b, 0.0, { ptIDs[0], ptIDs[1], ptIDs[2] });
    if (ptIDs[3] < 0) {
    	return false;
    }

    return true;
}

void getSInnerAndOuter(const trimesh::point& q, float r, float err, const C_RTree* rtree, std::vector<int>& ret, int excl_id_list[] = {}, int excl_id_num = 0) {

    rtree->nn_sphere_range(q, sq(r + err), err * 2.0, ret, excl_id_list, excl_id_num);

    std::unordered_set<int> new_excl_id_set;
    for (int i = 0; i < excl_id_num; i++) {
        new_excl_id_set.insert(excl_id_list[i]);
    }
    for (auto &v: ret) {
        new_excl_id_set.insert(v);
    }

    rtree->range_sphere_min_max(q, r - err, r + err, ret, new_excl_id_set);
}


void quickDonutRange(const trimesh::point& m, const trimesh::point& n, float r_donut, float err, const trimesh::TriMesh* mesh, const C_RTree* rtree, bool denseFlag, std::vector<int>& ret)
{
    int qDonutID = quickDonutNN(m, n, r_donut, mesh, rtree, denseFlag);
    if (qDonutID < 0 || qDonutID >= mesh->vertices.size()) {
        return;
    }

    const trimesh::point& qDonutPt = mesh->vertices[qDonutID];
    float qDonutDist = sqrt(donutDist2(qDonutPt, m, n, r_donut));
    float L = dist(qDonutPt, m);

    std::vector<int> donut_range_cand;
    rtree->range_sphere_min_max(m, std::max(L - err, 0.0f), L + err, donut_range_cand);

    for (auto &i: donut_range_cand) {
        const trimesh::point& i_pt = mesh->vertices[i];
        float donut_dist = sqrt(donutDist2(i_pt, m, n, r_donut));
        if (donut_dist <= qDonutDist + err * 2.0) {
            ret.push_back(i);
        }
    }
}

int calDonutCandidates(int q1ID, float r, float epsilon, const trimesh::TriMesh* mesh, const C_RTree* rtree, bool innerOnly, bool denseFlag, bool verbose, std::vector<int>& retList)
{
    int retSize = 0;
    float err = epsilon * 2.0;

    std::vector<int> rangeQ2, rangeQ3, rangeQ4;

    const trimesh::point& q1 = mesh->vertices[q1ID];

    if (innerOnly) {
        rtree->range_sphere_min_max(q1, r - err, r + err, rangeQ2);
    } else {
        getSInnerAndOuter(q1, r, err, rtree, rangeQ2);
    }

    for (auto &iQ2: rangeQ2) {

        if (iQ2 < 0)
            continue;

        // std::cout << iQ2 << std::endl;

        const trimesh::point& q2 = mesh->vertices[iQ2];

        auto m = middle_pt(q1, q2);
        auto n = q1 - q2;

        const float r_donut = dist(q1, m) * SQRT_3;
        quickDonutRange(m, n, r_donut, err, mesh, rtree, denseFlag, rangeQ3);

        for (auto &iQ3: rangeQ3) {

            if (iQ3 < 0)
                continue;

            // std::cout << iQ2 << "->" << iQ3 << std::endl;

            const trimesh::point& q3 = mesh->vertices[iQ3];

            auto a = getNNPtOnDonut(q3, m, n, r_donut);
            auto b = getPtBOnDonut(a, m, n);

            int excl_q123[3] = { q1ID, iQ2, iQ3 };
            rtree->nn_sphere_range(b, 0.0, err * 2.0, rangeQ4, excl_q123, 3);

            for (auto &iQ4: rangeQ4) {

                if (iQ4 < 0)
                    continue;

                // std::cout << iQ2 << "->" << iQ3 << "->" << iQ4 << std::endl;

                retList.push_back(iQ2);
                retList.push_back(iQ3);
                retList.push_back(iQ4);

                retSize++;
            }
            rangeQ4.clear();
        }
        rangeQ3.clear();
    }
    rangeQ2.clear();

    return retSize;
}

void calRd(const int ptIDs[4], const trimesh::TriMesh* mesh, float rd[RD_NUM])
{
	for (int i = 0; i < 4; i++) {
		assert(ptIDs[i] >= 0);
		assert(ptIDs[i] < mesh->vertices.size());
	}

	rd[0] = dist(mesh->vertices[ptIDs[0]], mesh->vertices[ptIDs[1]]);
	rd[1] = dist(mesh->vertices[ptIDs[0]], mesh->vertices[ptIDs[2]]);
	rd[2] = dist(mesh->vertices[ptIDs[0]], mesh->vertices[ptIDs[3]]);
	rd[3] = dist(mesh->vertices[ptIDs[1]], mesh->vertices[ptIDs[2]]);
	rd[4] = dist(mesh->vertices[ptIDs[1]], mesh->vertices[ptIDs[3]]);
	rd[5] = dist(mesh->vertices[ptIDs[2]], mesh->vertices[ptIDs[3]]);

	// might be different if e.g. color is used
}

class CFDetector
{
public:
    struct SingleCollTree {
        static const int RSTREE_SCALE = 1e5;
        static const int TREE_BUILD_THRESHOLD = 10;
        std::unordered_map<int, const trimesh::point*> insertedPtIds;
        CollTree* tree;
        std::unordered_map<int, bool> checkedSet;
        SingleCollTree() {
            this->tree = nullptr;
        }
        void insertTree(int pt_id, const trimesh::point* pt, bool verbose = false) {
            int pt_box_min[3], pt_box_max[3];
            for (int i = 0; i < 3; i++) {
                pt_box_min[i] = pt_box_max[i] = int((*pt)[i] * float(RSTREE_SCALE));
            }
            this->tree->Insert(pt_box_min, pt_box_max, pt_id);
            if (verbose) {
                std::cout << "Insert pt #" << pt_id << " into tree with box";
                for (int i = 0; i < 3; i++) {
                    std::cout << " [" << pt_box_min[i] << ", " << pt_box_max[i] << "]";
                }
                std::cout << std::endl;
            }
        }
        bool insert(int pt_id, const trimesh::point* pt, bool verbose = false) {
            if (this->insertedPtIds.find(pt_id) != this->insertedPtIds.end()) {
                return false;
            }

            this->insertedPtIds[pt_id] = pt;
            if (this->tree != nullptr) {
                this->insertTree(pt_id, pt, verbose);
            } else {
                if (this->insertedPtIds.size() > TREE_BUILD_THRESHOLD) {
                    if (verbose) {
                        std::cout << "Start a tree" << std::endl;
                    }
                    this->tree = new CollTree;
                    for (auto &i: this->insertedPtIds) {
                        insertTree(i.first, i.second, verbose);
                    }
                }
            }
            return true;
        }
        bool checkColl(int q_id, const trimesh::point* q, float d, float d2, bool& cacheRead, bool verbose = false) {
            if (checkedSet.find(q_id) != checkedSet.end()) {
                cacheRead = true;
                // if (verbose) {
                //     std::cout << "Read from cache: " << (checkedSet[q_id] ? "left" : "filtered") << std::endl;
                // }
                return checkedSet[q_id];
            }

            bool ret = false;
            if (this->tree != nullptr) {
                int q_box_min[3], q_box_max[3];
                for (int i = 0; i < 3; i++) {
                    q_box_min[i] = int(((*q)[i] - d) * float(RSTREE_SCALE));
                    q_box_max[i] = int(((*q)[i] + d) * float(RSTREE_SCALE));
                }
                ret = this->tree->YoNSearch(q_box_min, q_box_max); // TODO: use a fast SphereSearch instead, exclude itself (no need)?
                if (verbose) {
                    std::cout << "Search with box";
                    for (int i = 0; i < 3; i++) {
                        std::cout << " [" << q_box_min[i] << ", " << q_box_max[i] << "]";
                    }
                    std::cout << ": " << (ret ? "left" : "filtered") << std::endl;
                }
            } else {
                for (auto &p: this->insertedPtIds) {
                    // cout << " check dist=" << sq_dist(v, q);
                    if (dist2(*(p.second), *q) <= d2) {
                        ret = true;
                        break;
                    }
                }
                if (verbose) {
                    std::cout << "Search linearly: " << (ret ? "left" : "filtered") << std::endl;
                }
            }
            checkedSet[q_id] = ret;
            return ret;
        }
        void printInsertedPts(std::ostream &os) {
            for (auto &p: insertedPtIds) {
                os << p.first << *(p.second) << " ";
            }
        }
        void finish() {
            if (this->tree != nullptr) {
                this->tree->SortDim0();
            }
        }
        virtual ~SingleCollTree() {
            if (this->tree != nullptr) {
                delete this->tree;
            }
        }
    };

    std::unordered_map<int, SingleCollTree*> collMap;
    float nnDist;
    float nnDist2;

    int numCacheRead;

    CFDetector(float _nnDist) {
        this->nnDist = _nnDist;
        this->nnDist2 = _nnDist * _nnDist;

        numCacheRead = 0;
    }

    bool insert(int db_id, int pt_id, const trimesh::point* pt, bool verbose = false) {
        if (collMap.find(db_id) == collMap.end()) {
            collMap[db_id] = new SingleCollTree;
        }
        if (verbose) {
            std::cout << "Insert for DB #" << db_id << ", pt #" << pt_id << " " << *pt << std::endl;
        }

        return collMap[db_id]->insert(pt_id, pt, verbose);
    }

    bool checkColl(int db_id, int q_id, const trimesh::point* q, bool verbose = false) {
        if (verbose) {
            std::cout << "Check for DB #" << db_id << ", pt #" << q_id << " " << *q << ": ";
        }

        if (collMap.find(db_id) == collMap.end()) {
            if (verbose) {
                std::cout << "DB not exists, thus filtered" << std::endl;
            }
            return false;
        }

        bool cacheRead = false;
        bool ret = collMap[db_id]->checkColl(q_id, q, this->nnDist, this->nnDist2, cacheRead, verbose);
        if (cacheRead) {
            this->numCacheRead++;
        }
        return ret;
    }

    void finish() {
        for (auto &p: collMap) {
            p.second->finish();
        }
    }

    virtual ~CFDetector() {
        for (auto &p: collMap) {
            delete p.second;
        }
    }
};

#endif