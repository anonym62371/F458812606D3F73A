#ifndef __SUPER4PCS_H_
#define __SUPER4PCS_H_

#include <iostream>
#include <utility>

#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "TriMesh.h"
#include "KDtree.h"
// #include "ICP.h"
#include "TmpICP.hpp"

#include "MeshSet.hpp"
#include "UtilTime.hpp"

struct Params4PCS
{
	float delta = 0;
	bool prob = false;
	float epsilon = 0;
	int numIterations = 1;
	bool denseFlag = false;
	bool stopImmed = false;
	float fixedEpsilonDelta = -1;
};

void printParams4PCS(const Params4PCS& params)
{
	std::cout << "Parameters:" << std::endl;
	std::cout << "delta: " << params.delta << std::endl;
	std::cout << "prob: " << params.prob << std::endl;
	std::cout << "epsilon: " << params.epsilon << std::endl;
}

struct RetMatch4PCS
{
	uint32_t id;
	trimesh::xform xf;
	float dist;
	RetMatch4PCS(uint32_t _id, const trimesh::xform& _xf, float _dist) {
		this->id = _id;
		this->xf = _xf;
		this->dist = _dist;
	}
};

typedef std::vector<RetMatch4PCS> RetMatches4PCS;

struct RetStat4PCS
{
	int numItr = 0;

	float veriTime = 0;
	float veriCount = 0;
	float queryTime = 0;

	float propTime = 0;
};

const float COPLANAR_THRESHOLD = 0.01;
const float MATCH4PCS_ANGLE_THRESHOLD = 0.01;

class Super4PCS
{
private:
	MeshSet* dbMeshes;
	MeshSet* qMeshes;
	PairDist* pd;

	Params4PCS params4PCS;

	struct CongrBase {
	    int bp[4];
	    int db_id;
	    CongrBase(int bp1, int bp2, int bp3, int bp4, int db_id) {
	        this->bp[0] = bp1;
	        this->bp[1] = bp2;
	        this->bp[2] = bp3;
	        this->bp[3] = bp4;
	        this->db_id = db_id;
	    }
	};

    struct ExtractedPair {
        int p1, p2;
        trimesh::point e;
        ExtractedPair(int i1, int i2, const trimesh::TriMesh* mesh, float r) {
            this->p1 = i1;
            this->p2 = i2;
            auto& pt1 = mesh->vertices[i1];
            auto& pt2 = mesh->vertices[i2];
            e = (pt2 - pt1) * r + pt1;
        }
    };

	bool cal_intersection(int ran1, int ran2, int ran3, int ran4, const trimesh::TriMesh* qMesh, float& la, float& mu, float& angle);
	bool test_coplanar(int* ran, const trimesh::TriMesh* qMesh, float& d12_q, float& d34_q, float& la, float& mu, float& angle);

    int match_4pcs(const trimesh::TriMesh* mesh_p, const C_RTree* r_tree, std::vector<ExtractedPair>& ep12, std::unordered_map<int, std::vector<ExtractedPair*>>& ep34_map,
        float d34_q, float mu, float angle, int db_id, std::vector<CongrBase>& ret);

    int match_base_in_db_mesh_dp(float d12_q, float d34_q, float la, float mu, float angle, std::vector<CongrBase>& ret);
    
    void get_ep(const trimesh::TriMesh* mesh_p, const C_RTree* r_tree, float d, float ratio, std::vector<ExtractedPair>& ep_list);

    int match_base_in_db_mesh(const trimesh::TriMesh* mesh_p, const C_RTree* r_tree, float d12_q, float d34_q, float la, float mu, float angle, int db_id, std::vector<CongrBase>& ret);

	int oneIteration(const trimesh::TriMesh* qMesh, const C_RTree* qTree, RetMatches4PCS& retMatches4PCS, RetStat4PCS* retStat4PCS = nullptr, bool verbose = false);

public:
	Super4PCS(MeshSet* _dbMeshes, MeshSet* _qMeshes, PairDist* _pd = nullptr) {
		this->dbMeshes = _dbMeshes;
		this->qMeshes = _qMeshes;
		this->pd = _pd;
	}

	virtual ~Super4PCS() {
	}

	int execute(const Params4PCS& params, RetMatches4PCS& retMatches4PCS, RetStat4PCS* retStat4PCS = nullptr, bool verbose = false);
};

bool Super4PCS::cal_intersection(int ran1, int ran2, int ran3, int ran4, const trimesh::TriMesh* qMesh, float& la, float& mu, float& angle)
{
    auto& q1 = qMesh->vertices[ran1];
    auto& q2 = qMesh->vertices[ran2];
    auto& q3 = qMesh->vertices[ran3];
    auto& q4 = qMesh->vertices[ran4];
    mu = ((q3.x - q1.x) * (q2.y - q1.y) - (q3.y - q1.y) * (q2.x - q1.x))
       / ((q4.y - q3.y) * (q2.x - q1.x) - (q4.x - q3.x) * (q2.y - q1.y));

    // cout << mu << endl;
    if (mu < 0 || mu > 1) {
        return false;
    }

    auto e = q3 + (q4 - q3) * mu;

    la = (e.x - q1.x) / (q2.x - q1.x);

    // cout << la << endl;
    if (la < 0 || la > 1) {
        return false;
    }

    auto eq1 = q1 - e;
    auto eq3 = q3 - e;
    angle = (trimesh::dot(eq1, eq3)) / (trimesh::len(eq1) * trimesh::len(eq3));

    return true;
}

bool Super4PCS::test_coplanar(int* ran, const trimesh::TriMesh* qMesh, float& d12_q, float& d34_q, float& la, float& mu, float& angle)
{
    const trimesh::point* bq[4];
    int m = qMesh->vertices.size();
    int t = 0;
    bool found_coplanar = false;

    const int MAX_TRIAL = 50000;

    while (t < MAX_TRIAL && !found_coplanar) {
        // gen 3 different random indices in [0, m)
        ran[0] = rand() % m; ran[1] = 0; ran[2] = 0;
        do {
            ran[1] = rand() % m;
        } while (ran[1] == ran[0]);
        do {
            ran[2] = rand() % m;
        } while (ran[2] == ran[0] || ran[2] == ran[1]);

        // ran[0] = 430;
        // ran[1] = 191;
        // ran[2] = 322;
        // cout << "super4pcs rand points reset to " << ran[0] << ", " << ran[1] << ", " << ran[2] << endl;

        // get the first 3 base points
        for (int i = 0; i < 3; i++) {
            bq[i] = &(qMesh->vertices[ran[i]]);
        }
        
        // calculate the equation plane
        float a1 = bq[1]->x - bq[0]->x;
        float b1 = bq[1]->y - bq[0]->y;
        float c1 = bq[1]->z - bq[0]->z;
        float a2 = bq[2]->x - bq[0]->x;
        float b2 = bq[2]->y - bq[0]->y;
        float c2 = bq[2]->z - bq[0]->z;
        float a = b1 * c2 - b2 * c1;
        float b = a2 * c1 - a1 * c2;
        float c = a1 * b2 - b1 * a2;
        float d = (- a * bq[0]->x - b * bq[0]->y - c * bq[0]->z);

        float div = sqrt(a * a + b * b + c * c);

        // look for the fourth point
        for (int i = 0; i < m; i++) {
            if (i == ran[0] || i == ran[1] || i == ran[2]) // skip the repeats
                continue;

            bq[3] = &(qMesh->vertices[i]);
            float dist = (a * bq[3]->x + b * bq[3]->y + c * bq[3]->z + d) / div;

            if((abs(dist) <= COPLANAR_THRESHOLD) && this->cal_intersection(ran[0], ran[1], ran[2], i, qMesh, la, mu, angle)) {
                found_coplanar = true;
                ran[3] = i;
                break;
            }
        }

        if (found_coplanar) {
        	break;
        }

        t++;
    }

    if (!found_coplanar) {
        std::cout << "coplanar not found" << std::endl;
        return false;
    }

    d12_q = trimesh::dist(*(bq[0]), *(bq[1]));
    d34_q = trimesh::dist(*(bq[2]), *(bq[3]));

    std::cout << "coplanar found: " << ran[0] << " " << ran[1] << " " << ran[2] << " " << ran[3] << std::endl;
    std::cout << "d12: " << d12_q << std::endl;
    std::cout << "d34: " << d34_q << std::endl;
    std::cout << "la: " << la << std::endl;
    std::cout << "mu: " << mu << std::endl;
    std::cout << "angle: " << angle << std::endl;

    return true;
}

int Super4PCS::match_4pcs(const trimesh::TriMesh* mesh_p, const C_RTree* r_tree, std::vector<ExtractedPair>& ep12, std::unordered_map<int, std::vector<ExtractedPair*>>& ep34_map,
    float d34_q, float mu, float angle, int db_id, std::vector<CongrBase>& ret)
{
    float err = std::max(2.0 * this->params4PCS.epsilon, 0.001);

    float l34 = d34_q * mu;

    std::vector<int> intersection_ret;
    int num_matched = 0;

    for (auto &s12: ep12) {

        auto e12 = s12.e;
        auto& p1 = mesh_p->vertices[s12.p1];

        intersection_ret.clear();
        r_tree->range_sphere_dist_err(e12, l34, err, intersection_ret);

        for (auto &i: intersection_ret) {

            auto& iPt = mesh_p->vertices[i];
            auto e12_iPt = iPt - e12;
            auto e12_p1 = p1 - e12;
            double cosIPt = (dot(e12_iPt, e12_p1)) / (len(e12_iPt) * len(e12_p1));

            if (abs(cosIPt - angle) > MATCH4PCS_ANGLE_THRESHOLD) {
                continue;
            }

            auto got = ep34_map.find(i);
            if (got == ep34_map.end()) {
                // not found
                continue;
            }

            for (auto &s34: got->second) {
                float dist = trimesh::dist(e12, s34->e);

                if (dist <= err) {
                    ret.emplace_back(s12.p1, s12.p2, s34->p1, s34->p2, db_id);
                    num_matched++;
                }
            }
        }
    }

    return num_matched;
}

int Super4PCS::match_base_in_db_mesh_dp(float d12_q, float d34_q, float la, float mu, float angle, std::vector<CongrBase>& ret)
{
    float err = std::max(2.0 * this->params4PCS.epsilon, 0.001);

    std::vector<PairDistEntry*> pd_ret12;
    this->pd->searchEntries(std::max(0.0f, d12_q - err), d12_q + err, pd_ret12);
    std::vector<PairDistEntry*> pd_ret34;
    this->pd->searchEntries(std::max(0.0f, d34_q - err), d34_q + err, pd_ret34);

    std::cout << "Total number of matched pairs-12: " << pd_ret12.size() << std::endl;
    std::cout << "Total number of matched pairs-34: " << pd_ret34.size() << std::endl;

    // for (auto &v: pd_ret12) {
    //     cout << v->m_ep1 << " " << v->m_ep2 << " " << v->m_db_id << " " << v->m_dist << endl;
    // }
    // for (auto &v: pd_ret34) {
    //     cout << v->m_ep1 << " " << v->m_ep2 << " " << v->m_db_id << " " << v->m_dist << endl;
    // }

    std::unordered_map<int, std::vector<ExtractedPair>> ep12_map;
    std::unordered_map<int, std::vector<ExtractedPair>> ep34_map;
    for (auto &pde: pd_ret12) {
        if (ep12_map.find(pde->m_db_id) == ep12_map.end()) {
            std::vector<ExtractedPair> new_eps;
            ep12_map[pde->m_db_id] = new_eps;
        }
        ep12_map[pde->m_db_id].emplace_back(pde->m_ep1, pde->m_ep2, dbMeshes->getMesh(pde->m_db_id).mesh /*qc.db_meshes.get_mesh(pde->m_db_id)*/, la);
        ep12_map[pde->m_db_id].emplace_back(pde->m_ep2, pde->m_ep1, dbMeshes->getMesh(pde->m_db_id).mesh /*qc.db_meshes.get_mesh(pde->m_db_id)*/, la);
    }
    for (auto &pde: pd_ret34) {
        if (ep34_map.find(pde->m_db_id) == ep34_map.end()) {
            std::vector<ExtractedPair> new_eps;
            ep34_map[pde->m_db_id] = new_eps;
        }
        ep34_map[pde->m_db_id].emplace_back(pde->m_ep1, pde->m_ep2, dbMeshes->getMesh(pde->m_db_id).mesh /* qc.db_meshes.get_mesh(pde->m_db_id)*/, mu);
        ep34_map[pde->m_db_id].emplace_back(pde->m_ep2, pde->m_ep1, dbMeshes->getMesh(pde->m_db_id).mesh /* qc.db_meshes.get_mesh(pde->m_db_id)*/, mu);
    }

    int num_matched;

    for (auto &v: ep12_map) {
        int db_id = v.first;
        if (ep34_map.find(db_id) == ep34_map.end()) {
            continue;
        }
        auto ep12 = v.second;
        auto& ep34 = ep34_map[db_id];

        auto mesh_p = dbMeshes->getMesh(db_id).mesh;
        auto r_tree = dbMeshes->getMesh(db_id).rtree;

        std::unordered_map<int, std::vector<ExtractedPair*>> ep34_by_p1_map;
        for (int i = 0; i < ep34.size(); i++) {
            if (ep34_by_p1_map.find(ep34[i].p1) == ep34_by_p1_map.end()) {
                std::vector<ExtractedPair*> new_list;
                ep34_by_p1_map[ep34[i].p1] = new_list;
            }
            ep34_by_p1_map[ep34[i].p1].push_back(&ep34[i]);
        }

        // cout << "Process DB #" << db_id << endl;
        int num_matched_ind = match_4pcs(mesh_p, r_tree, ep12, ep34_by_p1_map, d34_q, mu, angle, db_id, ret);
        // cout << "Matched with: " << num_matched_ind << endl;

        num_matched += num_matched_ind;
    }

    return num_matched;
}

void Super4PCS::get_ep(const trimesh::TriMesh* mesh_p, const C_RTree* r_tree, float d, float ratio, std::vector<ExtractedPair>& ep_list)
{
    int n = mesh_p->vertices.size();
    // double err = max(0.01, 2 * qc.epsilon); // why this is 0.01
    float err = std::max(2.0 * this->params4PCS.epsilon, 0.01);
    std::vector<int> intersection_ret;

    for (int i = 0; i < n; i++) {

        auto& p = mesh_p->vertices[i];

        // #ifdef _CLR
        //     auto p_color = get_p_color(p, mesh_p);
        //     if (!(same_color(&p_color, qc.bq_clr[0]) || same_color(&p_color, qc.bq_clr[1]) || same_color(&p_color, qc.bq_clr[2]) || same_color(&p_color, qc.bq_clr[3]))) {
        //         continue;
        //     }
        // #endif

        intersection_ret.clear();
        r_tree->range_sphere_dist_err(p, d, err, intersection_ret);

        for (auto &r: intersection_ret) {
            ep_list.emplace_back(i, r, mesh_p, ratio);
        }
    }
}

int Super4PCS::match_base_in_db_mesh(const trimesh::TriMesh* mesh_p, const C_RTree* r_tree, float d12_q, float d34_q, float la, float mu, float angle, int db_id, std::vector<CongrBase>& ret)
{
    float err = std::max(2.0 * this->params4PCS.epsilon, 0.001);

    std::vector<ExtractedPair> ep12;
    std::vector<ExtractedPair> ep34;
    get_ep(mesh_p, r_tree, d12_q, la, ep12);
    get_ep(mesh_p, r_tree, d34_q, mu, ep34);

    // std::cout << "Num of matching 12: " << ep12.size() << std::endl;
    // std::cout << "Num of matching 34: " << ep34.size() << std::endl;

    std::unordered_map<int, std::vector<ExtractedPair*>> ep34_map;
    for (int i = 0; i < ep34.size(); i++) {
        if (ep34_map.find(ep34[i].p1) == ep34_map.end()) {
            std::vector<ExtractedPair*> new_list;
            ep34_map[ep34[i].p1] = new_list;
        }
        ep34_map[ep34[i].p1].push_back(&ep34[i]);
    }

    return match_4pcs(mesh_p, r_tree, ep12, ep34_map, d34_q, mu, angle, db_id, ret);
}

int Super4PCS::oneIteration(const trimesh::TriMesh* qMesh, const C_RTree* qTree, RetMatches4PCS& retMatches4PCS, RetStat4PCS* retStat4PCS, bool verbose)
{
    if (retStat4PCS) timer_start();

    int ran[4];
    float la, mu;
    float angle; // cos(theta)
    float d12_q, d34_q; // for super4pcs
    float l, d;
    float d1, d2; // for donut-basic

    bool found_coplanar = this->test_coplanar(ran, qMesh, d12_q, d34_q, la, mu, angle);
    /* other variants */
    // if (qc.algorithm == "donut-basic") {
    //     found_coplanar = donut_select(ran, d1, d2, l, la, mu, angle);
    // } else if (qc.algorithm == "superg4pcs-adv" || qc.algorithm == "superg4pcs-basic") {
    //     found_coplanar = test_3d_coplanar(ran, d, l, la, mu, angle);
    // } else {
    //     found_coplanar = test_coplanar(ran, d12_q, d34_q, la, mu, angle);
    // }

    if (!found_coplanar) {
        if (retStat4PCS) retStat4PCS->queryTime += timer_end();
        return 0;
    }

    std::cout << "Select rand query points: " << ran[0] << " " << ran[1] << " " << ran[2] << " " << ran[3] << std::endl;

    std::vector<CongrBase> retCongrBase;

    if (this->pd) {
        match_base_in_db_mesh_dp(d12_q, d34_q, la, mu, angle, retCongrBase);
    } else {
        for (int i = 0; i < dbMeshes->size(); i++) {
            // cout << "Match base in db #" << i << endl;
            int ret_size = match_base_in_db_mesh(dbMeshes->getMesh(i).mesh, dbMeshes->getMesh(i).rtree, d12_q, d34_q, la, mu, angle, i, retCongrBase);
            // cout << "Returned matched: " << ret_size << endl;
        }
    }

    std::cout << "Number of returned congr-bases: " << retCongrBase.size() << std::endl;

    // for (auto &cb: retCongrBase) {
    //     std::cout << "Matched with DB #" << cb.db_id << ": " << cb.bp[0] << " " << cb.bp[1] << " " << cb.bp[2] << " " << cb.bp[3] << std::endl;
    // }

    /* other variants */
    // if (qc.algorithm == "super4pcs-adv") {
    //     match_base_in_db_mesh_dp(d12_q, d34_q, la, mu, angle, retCongrBase);
    // } else if (qc.algorithm == "superg4pcs-adv") {
    //     match_base_in_db_mesh_gen_dp(d, l, la, mu, angle, retCongrBase);
    // } else {
    //     for (int i = 0; i < qc.db_meshes.size(); i++) {
    //         // cout << "Match base in db#" << i << endl;
    //         int ret_size;
    //         if (qc.algorithm == "donut-basic") {
    //             ret_size = match_base_in_db_mesh_donut(qc.db_meshes.get_mesh(i), qc.db_rtrees[i], d1, d2, l, la, mu, angle, i, retCongrBase);
    //         } else if (qc.algorithm == "superg4pcs-basic") {
    //             ret_size = match_base_in_db_mesh_gen(qc.db_meshes.get_mesh(i), qc.db_rtrees[i], d, l, la, mu, angle, i, retCongrBase);
    //         } else if (qc.algorithm == "super4pcs-basic") {
    //             ret_size = match_base_in_db_mesh(qc.db_meshes.get_mesh(i), qc.db_rtrees[i], d12_q, d34_q, la, mu, angle, i, retCongrBase);
    //         }
    //         // cout << "Returned matched: " << ret_size << endl;
    //     }
    // }

    const trimesh::point* qPts[4];
    for (int i = 0; i < 4; i++) {
        qPts[i] = &(qMesh->vertices[ran[i]]);
    }

    /* perform query */
    const trimesh::point* pPts[4];
    trimesh::xform xf;
    int numMatched = 0;
    int numVerified = 0;

    std::unordered_set<int> verifiedIDSet;

    for (auto &cb: retCongrBase) {

        if (verifiedIDSet.find(cb.db_id) != verifiedIDSet.end()) {
            // std::cout << "Filtered: " << cb.db_id << std::endl;
            continue;
        }

        for (int i = 0; i < 4; i++) {
            pPts[i] = &(this->dbMeshes->getMesh(cb.db_id).mesh->vertices[cb.bp[i]]);
        }

        if (retStat4PCS) timer_start();

        cal_trans(qPts, pPts, 4, xf);

        tmp_icp_v2v(qMesh, this->dbMeshes->getMesh(cb.db_id).mesh, xf, this->dbMeshes->getMesh(cb.db_id).kdtree);

        float pairDist2 = this->dbMeshes->calDist2(qMesh, cb.db_id, &xf, sq(this->params4PCS.delta) * float(qMesh->vertices.size()));
        // std::cout << "Returned dist: " << pairDist2 << std::endl;

        numVerified ++;

        if (retStat4PCS) {
            retStat4PCS->veriTime += timer_end();
            retStat4PCS->veriCount ++;
        }

        if (pairDist2 >= 0.0f) {
            float convertedDist = sqrt(pairDist2 / float(qMesh->vertices.size())) / qMesh->bsphere.r / 2.0;

            std::cout << "Matched with database mesh object " << this->dbMeshes->getMesh(cb.db_id).filename;
            std::cout << " of distance " << convertedDist << std::endl;
            std::cout << xf;
            // std::cout << "ICP returned dist: " << pairDistICP << std::endl;

            retMatches4PCS.emplace_back(cb.db_id, xf, convertedDist);
            numMatched++;

            verifiedIDSet.insert(cb.db_id);

            if (this->params4PCS.stopImmed) {
                if (verbose) std::cout << "Return immediately" << std::endl;
                break;
            }
        }

    }

    if (retStat4PCS) retStat4PCS->queryTime += timer_end();

    std::cout << "Number of total pairs verified: " << numVerified << std::endl;

    return numMatched;
}


int Super4PCS::execute(const Params4PCS& params, RetMatches4PCS& retMatches4PCS, RetStat4PCS* retStat4PCS, bool verbose)
{
	this->params4PCS = params;

	int numItr = 1;
	if (this->params4PCS.prob) {
		numItr = 8;
	}

	auto qMesh = qMeshes->getMesh(0).mesh;
	auto qRtree = qMeshes->getMesh(0).rtree;
	// auto qKd = qMeshes->getMesh(0).kdtree;

	if (!this->params4PCS.prob) {
		if (this->params4PCS.fixedEpsilonDelta > 0) {
			this->params4PCS.epsilon = this->params4PCS.fixedEpsilonDelta * sqrt(float(qMesh->vertices.size()));
		} else {
			this->params4PCS.epsilon = this->params4PCS.delta * sqrt(float(qMesh->vertices.size()));
		}
	}
	std::cout << "Epsilon: " << this->params4PCS.epsilon << std::endl;

	// printParams4PCS(this->params4PCS);

	int matchedSize = 0;
	for (int i = 0; i < numItr; i++) {
		if (retStat4PCS) retStat4PCS->numItr++;

		matchedSize = this->oneIteration(qMesh, qRtree, retMatches4PCS, retStat4PCS, verbose);
		if (matchedSize > 0) {
			break;
		}
	}

	// derived stat:
	if (retStat4PCS) {
		retStat4PCS->propTime = retStat4PCS->queryTime - retStat4PCS->veriTime;
	}

	return matchedSize;
}


#endif