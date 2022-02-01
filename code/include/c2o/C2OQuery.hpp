#ifndef __C2OQUERY_H_
#define __C2OQUERY_H_

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
#include "C2OIndex.hpp"
#include "UtilTime.hpp"
#include "UtilStats.hpp"

struct QueryParams
{
	float delta = 0;
	int hCF = 1;
	bool prob = false;
	float epsilon = 0;
	int numIterations = 1;
	bool innerOnly = false;
	bool denseFlag = false;
	bool stopImmed = false;
	float fixedEpsilonDelta = -1;
};

void printQueryParams(const QueryParams& params)
{
	std::cout << "Parameters:" << std::endl;
	std::cout << "delta: " << params.delta << std::endl;
	std::cout << "hCF: " << params.hCF << std::endl;
	std::cout << "prob: " << params.prob << std::endl;
	std::cout << "epsilon: " << params.epsilon << std::endl;
}

struct RetMatch
{
	uint32_t id;
	trimesh::xform xf;
	float dist;
	RetMatch(uint32_t _id, const trimesh::xform& _xf, float _dist) {
		this->id = _id;
		this->xf = _xf;
		this->dist = _dist;
	}
};

typedef std::vector<RetMatch> RetMatches;

struct RetStat
{
	int numItr = 0;

	float calDonutCandTimeWoCF = 0;
	float donutCandCountAvg = 0;

	float CFBuildTime = 0;
	
	float searchEntriesTime = 0;
	int searchEntriesCount = 0;
	float searchEntriesRetSizeAvg = 0;

	float pruningPower1Avg = 0;
	float pruningPower2Avg = 0;
	
	float CFCheckTime = 0;
	int CFCheckCount = 0;
	
	float veriTime = 0;
	float veriCount = 0;
	float queryTime = 0;

	float propTime = 0;
};


class C2OQuery
{
private:
	MeshSet* dbMeshes;
	MeshSet* qMeshes;
	C2OIndex* c2oIndex;

	QueryParams queryParams;

	const C2OVariantEnum c2oVariant = C2OV_DONUT;

	struct DonutRetList {
		int randID;
		std::vector<int> donutRetList;
		int donutRetSize;
	};

	int getDrl(const trimesh::TriMesh* qMesh, const C_RTree* qTree, int randID, bool verbose, DonutRetList& drl);

	void getRandIDSet(const trimesh::TriMesh* qMesh, std::unordered_set<int>& randIDSet);

	int oneIteration(const trimesh::TriMesh* qMesh, const C_RTree* qTree, const trimesh::KDtree* qKd, RetMatches& retMatches, RetStat* retStat = nullptr, bool verbose = false);

public:
	C2OQuery(MeshSet* _dbMeshes, MeshSet* _qMeshes, C2OIndex* _c2oIndex) {
		this->dbMeshes = _dbMeshes;
		this->qMeshes = _qMeshes;
		this->c2oIndex = _c2oIndex;
	}

	virtual ~C2OQuery() {
	}

	int execute(const QueryParams& params, RetMatches& retMatches, RetStat* retStat = nullptr, bool verbose = false);

};

int C2OQuery::getDrl(const trimesh::TriMesh* qMesh, const C_RTree* qTree, int randID, bool verbose, DonutRetList& drl)
{
	int ret = calDonutCandidates(randID,
		this->c2oIndex->m_indexParams.r,
		this->queryParams.epsilon,
		qMesh,
		qTree,
		this->queryParams.innerOnly,
		this->queryParams.denseFlag,
		verbose,
		drl.donutRetList
	);
	drl.donutRetSize = ret;
	drl.randID = randID;
	return ret;
}

void C2OQuery::getRandIDSet(const trimesh::TriMesh* qMesh, std::unordered_set<int>& randIDSet)
{
	const int MAX_RAND_TRIAL = 5;
	// const int MAX_RAND_TRIAL = 1;

	int counter = 0;
	// first try to select points near to the center
	while (randIDSet.size() < MAX_RAND_TRIAL && counter < MAX_RAND_TRIAL * 5) {
		counter++;
		
		int randID = rand() % qMesh->vertices.size();
		if (randIDSet.find(randID) != randIDSet.end()) {
			continue;
		}

		const trimesh::point& randPt = qMesh->vertices[randID];
		if (dist(randPt, qMesh->bsphere.center) > 0.8 * (qMesh->bsphere.r - this->c2oIndex->m_indexParams.r)) {
			continue;
		}

		randIDSet.insert(randID);
	}

	if (!randIDSet.empty()) {
		return;
	}

	// if the previous selection failed, we pick points randomly
	// assume there are at least MAX_RAND_TRIAL points in qMesh
	while (randIDSet.size() < MAX_RAND_TRIAL) {
		int randID = rand() % qMesh->vertices.size();
		if (randIDSet.find(randID) != randIDSet.end()) {
			continue;
		}
		randIDSet.insert(randID);
	}
}

int C2OQuery::oneIteration(const trimesh::TriMesh* qMesh, const C_RTree* qTree, const trimesh::KDtree* qKd, RetMatches& retMatches, RetStat* retStat, bool verbose)
{
	retMatches.clear();

	if (retStat) timer_start();
	if (retStat) timer_start();

	/* 1a). select a set of random points in query near center */
	std::unordered_set<int> randIDSet;
	getRandIDSet(qMesh, randIDSet);

	/* 1b). for each random point, find the donut list and sort by return list size */
	std::vector<DonutRetList> drls = std::vector<DonutRetList>(randIDSet.size());
	int iDrl = 0;
	for (auto &randID: randIDSet) {
		if (verbose) std::cout << "Get DRL for pt #" << randID << std::endl;
		this->getDrl(qMesh, qTree, randID, verbose, drls[iDrl]);
		iDrl++;
	}
	std::sort(drls.begin(), drls.end(), [](DonutRetList const& drl1, DonutRetList const& drl2) {
		return drl1.donutRetSize < drl2.donutRetSize;
	});
	if (verbose) {
		std::cout << "Trial list:";
		for (auto &drl: drls) {
			std::cout << " " << drl.randID << "(" << drl.donutRetSize << ")";
		}
		std::cout << std::endl;
	}

	if (retStat) retStat->calDonutCandTimeWoCF += timer_end();

	/* 1c). find the smallest non-empty donut ret list */
	int startingDrlFlag = 0;
	for (int i = 0; i < drls.size(); i++) {
		if (drls[i].donutRetSize <= 0) {
			startingDrlFlag++;
		} else {
			break;
		}
	}
	if (startingDrlFlag >= drls.size()) {
		std::cerr << "Selected random points lead to empty donut lists" << std::endl;
		return 0;
	}
	const DonutRetList* ptrSelDrl = &(drls[startingDrlFlag]);

	// // debug: manually set
	// DonutRetList drlManual;
	// // this->getDrl(qMesh, qTree, 81, verbose, drlManual);
	// this->getDrl(qMesh, qTree, 97, verbose, drlManual);
	// ptrSelDrl = &drlManual;
	// std::cout << "Manually set DRL of pt #" << ptrSelDrl->randID << std::endl;
	// // debug: manually set

	// std::cout << "Select random query pt #" << ptrSelDrl->randID << std::endl;
	// std::cout << "Returned donut candidates of number " << ptrSelDrl->donutRetSize << std::endl;

	if (retStat) {
		retStat->donutCandCountAvg
			= retStat->donutCandCountAvg / float(retStat->numItr) * float(retStat->numItr - 1)
			+ float(ptrSelDrl->donutRetSize) / float(retStat->numItr);
	}

	/* 1.5. necessary structures for RD-pair retrieval*/
	float err[RD_NUM];
	for (int i = 0; i < RD_NUM; i++) {
		err[i] = this->queryParams.epsilon * 2.0;
	}
	float qRd[RD_NUM];
	int qCandIDs[4] = { ptrSelDrl->randID }; // first is randID
	const trimesh::point* qPts[4] = { &(qMesh->vertices[qCandIDs[0]]) };

	if (retStat) timer_start();

	/* 2a). find the k-nn point IDs for CF */
	int retCF[queryParams.hCF];
	if (queryParams.hCF > 0) {
		qTree->knn_sphere(qMesh->vertices[ptrSelDrl->randID], 0.01, queryParams.hCF, retCF);
	}
	if (verbose) {
		if (queryParams.hCF > 0) {
			std::cout << "CF pt ID list by dist:";
			for (int i = 0; i < queryParams.hCF; i++) {
				std::cout << " " << retCF[i] << "(" << dist(*(qPts[0]), qMesh->vertices[retCF[i]]) << ")";
			}
			std::cout << std::endl;
		}
	}

	/* 2b). calculate the CF detectors */
	CFDetector* cfds[queryParams.hCF];
	for (int i = 0; i < queryParams.hCF; i++) {
		if (verbose) std::cout << "Get DRL for CF pt #" << retCF[i] << std::endl;
		DonutRetList* ptrDrlCF = nullptr;
		DonutRetList drlCFBackup;
		for (int j = 0; j < drls.size(); j++) {
			if (drls[j].randID == retCF[i]) {
				ptrDrlCF = &(drls[j]);
				if (verbose) std::cout << "Find in cache at pos " << j << std::endl;
				break;
			}
		}
		if (!ptrDrlCF) {
			this->getDrl(qMesh, qTree, retCF[i], verbose, drlCFBackup);
			ptrDrlCF = &drlCFBackup;
			if (verbose) std::cout << "Not find in cache, rebuild with candidates size: " << ptrDrlCF->donutRetSize << std::endl;
		}

		int qCandIDsCF[4] = { ptrDrlCF->randID }; // first is CF
		const auto q1CF = qMesh->vertices[qCandIDsCF[0]];
		float cfDist = dist(*(qPts[0]), q1CF) + this->queryParams.epsilon * 2.0;

		cfds[i] = new CFDetector(cfDist);
		// std::cout << "CF #" << i << " initialized with pt #" << ptrDrlCF->randID << " dist " << cfDist << std::endl;

		std::vector<const IndexEntry*> pairRetListCF;
		int totalCF = 0;
		for (int j = 0; j < ptrDrlCF->donutRetSize; j++) {
			qCandIDsCF[1] = ptrDrlCF->donutRetList[j*3];
			qCandIDsCF[2] = ptrDrlCF->donutRetList[j*3+1];
			qCandIDsCF[3] = ptrDrlCF->donutRetList[j*3+2];

			calRd(qCandIDsCF, qMesh, qRd);
			this->c2oIndex->searchEntries(qRd, err, pairRetListCF);

			for (auto &pair: pairRetListCF) {
				auto p1 = &(this->dbMeshes->getMesh(pair->dbID).mesh->vertices[pair->ptIDs[0]]);

				bool vInsert = false;
				// if (pair->dbID == 978) {
				// 	vInsert = true;
				// }
				if (cfds[i]->insert(pair->dbID, pair->ptIDs[0], p1, vInsert)) {
					totalCF++;
				}
			}

			pairRetListCF.clear();
		}

		cfds[i]->finish();
		// std::cout << "CF #" << i << " built with number of total CF points: " << totalCF << std::endl;
	}

	if (retStat) retStat->CFBuildTime += timer_end();
	
	// const int testID = 978;
	// std::cout << "Test view in CF: DB #" << testID << ": ";
	// cfds[0]->collMap[testID]->printInsertedPts(std::cout);
	// std::cout << std::endl;
	// std::cout << cfds[0]->collMap[testID]->tree->Count() << std::endl;


	/* perform query */
	const trimesh::point* pPts[4];
	std::vector<const IndexEntry*> pairRetList;
	trimesh::xform xf;
	int numMatched = 0;

	int numTotal = 0;
	int numCFLeft = 0;

	std::unordered_set<int> verifiedIDSet;

	for (int i = 0; i < ptrSelDrl->donutRetSize; i++) {
		qCandIDs[1] = ptrSelDrl->donutRetList[i*3];
		qCandIDs[2] = ptrSelDrl->donutRetList[i*3+1];
		qCandIDs[3] = ptrSelDrl->donutRetList[i*3+2];

		calRd(qCandIDs, qMesh, qRd);

		pairRetList.clear();

		if (retStat) timer_start();
		int nodeAccessed = 0;
		int leafAccessed = 0;
		this->c2oIndex->searchEntries(qRd, err, pairRetList, &nodeAccessed, &leafAccessed);
		// std::cout << "node_accessed: " << node_accessed << ", leaf_accessed: " << leaf_accessed << ", ret: " << pairRetList.size() << std::endl;
		if (retStat) {
			retStat->searchEntriesTime += timer_end();
			retStat->searchEntriesCount ++;
			float searchEntriesRetSize = pairRetList.size();
			retStat->searchEntriesRetSizeAvg
				= retStat->searchEntriesRetSizeAvg / float(retStat->searchEntriesCount) * float(retStat->searchEntriesCount - 1)
				+ searchEntriesRetSize / float(retStat->searchEntriesCount);

			float pruningPower1 = 1.0 - float(nodeAccessed + leafAccessed) / this->dbMeshes->totalVertices();
			float pruningPower2 = 1.0 - float(leafAccessed) / this->dbMeshes->totalVertices();

			retStat->pruningPower1Avg
				= retStat->pruningPower1Avg / float(retStat->searchEntriesCount) * float(retStat->searchEntriesCount - 1)
				+ pruningPower1 / float(retStat->searchEntriesCount);

			retStat->pruningPower2Avg
				= retStat->pruningPower2Avg / float(retStat->searchEntriesCount) * float(retStat->searchEntriesCount - 1)
				+ pruningPower2 / float(retStat->searchEntriesCount);
		}

		qPts[1] = &(qMesh->vertices[qCandIDs[1]]);
		qPts[2] = &(qMesh->vertices[qCandIDs[2]]);
		qPts[3] = &(qMesh->vertices[qCandIDs[3]]);

		for (auto &pair: pairRetList) {
			if (verifiedIDSet.find(pair->dbID) != verifiedIDSet.end()) {
				// std::cout << "Filtered: " << pair->dbID << std::endl;
				continue;
			}

			pPts[0] = &(this->dbMeshes->getMesh(pair->dbID).mesh->vertices[pair->ptIDs[0]]);

			numTotal++;

			bool verboseCF = false;
			// if (pair->dbID == testID && pair->ptIDs[0] == 97) {
			// 	verboseCF = true;
			// }
			bool left = true;
			if (retStat) timer_start();
			for (int i = 0; i < queryParams.hCF; i++) {
				if (!cfds[i]->checkColl(pair->dbID, pair->ptIDs[0], pPts[0], verboseCF)) {
					left = false;
					break;
				}
			}
			if (retStat) {
				retStat->CFCheckTime += timer_end();
				retStat->CFCheckCount ++;
			}
			if (!left) {
				continue; // CF coll check failed
			}

			numCFLeft++;

			pPts[1] = &(this->dbMeshes->getMesh(pair->dbID).mesh->vertices[pair->ptIDs[1]]);
			pPts[2] = &(this->dbMeshes->getMesh(pair->dbID).mesh->vertices[pair->ptIDs[2]]);
			pPts[3] = &(this->dbMeshes->getMesh(pair->dbID).mesh->vertices[pair->ptIDs[3]]);

			// if (pair->dbID != 592) {
			// 	continue;
			// }

			// std::cout << "(" << qCandIDs[0] << " " << qCandIDs[1] << " " << qCandIDs[2] << " " << qCandIDs[3] << ") <--> ";
			// std::cout << "(" << pair->ptIDs[0] << " " << pair->ptIDs[1] << " " << pair->ptIDs[2] << " " << pair->ptIDs[3] << ")" << std::endl;

			if (retStat) timer_start();

			cal_trans(qPts, pPts, 4, xf);

			// float startPairDist2 = this->dbMeshes->calDist2(qMesh, pair->dbID, &xf);
			// std::cout << "Starting dist: " << (sqrt(startPairDist2 / float(qMesh->vertices.size())) / qMesh->bsphere.r / 2.0) << std::endl;

			tmp_icp_v2v(qMesh, this->dbMeshes->getMesh(pair->dbID).mesh, xf, this->dbMeshes->getMesh(pair->dbID).kdtree);

			// trimesh::TriMesh::set_verbose(1);
			// float pairDistICP = 
			// ICP(this->dbMeshes->getMesh(pair->dbID).mesh, qMesh, trimesh::xform(), xf, this->dbMeshes->getMesh(pair->dbID).kdtree, qKd, 1);
			// trimesh::TriMesh::set_verbose(0);

			// float refinedPairDist2 = this->dbMeshes->calDist2(qMesh, pair->dbID, &xf);
			// std::cout << "Refined dist: " << (sqrt(refinedPairDist2 / float(qMesh->vertices.size())) / qMesh->bsphere.r / 2.0) << std::endl;

			// std::cout << "Refined XF:" << std::endl;
			// std::cout << xf << std::endl;

			float pairDist2 = this->dbMeshes->calDist2(qMesh, pair->dbID, &xf, sq(this->queryParams.delta) * float(qMesh->vertices.size()));
			// std::cout << "Returned dist: " << pairDist2 << std::endl;

			if (retStat) {
				retStat->veriTime += timer_end();
				retStat->veriCount ++;
			}

			if (pairDist2 >= 0.0f) {
				float convertedDist = sqrt(pairDist2 / float(qMesh->vertices.size())) / qMesh->bsphere.r / 2.0;

				std::cout << "Matched with database object " << this->dbMeshes->getMesh(pair->dbID).filename;
				std::cout << " of distance " << convertedDist << std::endl;
				std::cout << xf;
				// std::cout << "ICP returned dist: " << pairDistICP << std::endl;

				retMatches.emplace_back(pair->dbID, xf, convertedDist);
				numMatched++;

				verifiedIDSet.insert(pair->dbID);
				// for (auto &v: verifiedIDSet) {
				// 	std::cout << " " << v;
				// }
				// std::cout << std::endl;

				if (this->queryParams.stopImmed) {
					if (verbose) std::cout << "Return immediately" << std::endl;
					break;
				}
			}
		}

		if (this->queryParams.stopImmed && numMatched > 0) {
			break;
		}
	}

	if (retStat) retStat->queryTime += timer_end();

	// std::cout << "Number of total pairs to verify: " << numTotal << std::endl;
	// std::cout << "Number of pairs left after CF: " << numCFLeft << std::endl;

	for (int i = 0; i < queryParams.hCF; i++) {
		if (cfds[i]) {
			delete cfds[i];
		}
	}

	return numMatched;
}


int C2OQuery::execute(const QueryParams& params, RetMatches& retMatches, RetStat* retStat, bool verbose)
{
	this->queryParams = params;

	int numItr = 1;
	if (this->queryParams.prob) {
		numItr = (this->queryParams.hCF > 0) ? 20 : 8;
	}

	auto qMesh = qMeshes->getMesh(0).mesh;
	auto qRtree = qMeshes->getMesh(0).rtree;
	auto qKd = qMeshes->getMesh(0).kdtree;

	if (!this->queryParams.prob) {
		if (this->queryParams.fixedEpsilonDelta > 0) {
			this->queryParams.epsilon = this->queryParams.fixedEpsilonDelta * sqrt(float(qMesh->vertices.size()));
		} else {
			this->queryParams.epsilon = this->queryParams.delta * sqrt(float(qMesh->vertices.size()));
		}
	}
	// std::cout << "Epsilon: " << this->queryParams.epsilon << std::endl;

	if (queryParams.prob) {
		this->queryParams.innerOnly = true;
	}

	// printQueryParams(this->queryParams);

	int matchedSize = 0;
	for (int i = 0; i < numItr; i++) {
		if (retStat) retStat->numItr++;

		matchedSize = this->oneIteration(qMesh, qRtree, qKd, retMatches, retStat, verbose);
		if (matchedSize > 0) {
			break;
		}
	}

	// derived stat:
	if (retStat) {
		retStat->propTime = retStat->queryTime - retStat->veriTime;
	}

	return matchedSize;
}


#endif