#ifndef __C2OINDEX_H_
#define __C2OINDEX_H_

#include <iostream>
#include <fstream>
#include <string>
#include <cstdint>
#include <filesystem>
#include <limits>

#include "ProgressBar.hpp"
#include "UtilVec.hpp"
#include "UtilTime.hpp"
#include "C2O.hpp"
#include "CPP_RTree.hpp"
#include "MeshSet.hpp"

namespace fs = std::filesystem;

struct IndexParams {
	float r;
	float r2, r3;
	float angleThresh;
	bool denseFlag = false;
};

void printIndexParams(const IndexParams& params)
{
	std::cout << "Parameters:" << std::endl;
	std::cout << "r: " << params.r << std::endl;
}

struct IndexEntry {
	// uint64_t globalID;
	int32_t dbID;
	int32_t ptIDs[4];
	IndexEntry(/*uint64_t parGlobalID, */int32_t parDbID, int32_t parPtIDs[4]) {
		// globalID = parGlobalID;
		dbID = parDbID;
		for (int i = 0; i < 4; i++) {
			ptIDs[i] = parPtIDs[i];
		}
	}
	IndexEntry() {
		dbID = -1;
		for (int i = 0; i < 4; i++) {
			ptIDs[i] = -1;
		}
	}
	#ifdef _IDX3
		float unindexedSides[3];
	#endif
	bool isValid() {
		if (dbID < 0) {
			return false;
		}
		for (auto &v: ptIDs) {
			if (v < 0) {
				return false;
			}
		}
		return true;
	}
};
// currently, to save space, all use 32-bit integers, and in the future, there might be 64-bit need
// in the future, support the index of other attributes

class C2OIndex
{
private:
	// struct TechFields {
	// 	uint32_t numEntries;
	// 	uint8_t enableIdx3;
	// 	uint8_t enableClr;
	// }; // For IO

	MeshSet* m_meshSet;

	std::vector<IndexEntry> m_entryList;

	IndexTree m_indexTree;

	const C2OVariantEnum c2oVariant = C2OV_DONUT;

	static const int RSTREE_SCALE = 1e5;

	// TODO: if rtree = null? then how: write some brute-force solutions
	bool calEntry(int ptIDs[4], const trimesh::TriMesh* mesh, const C_RTree* rtree, bool verbose = false) {
		switch(this->c2oVariant) {
			case C2OV_DONUT: return calDonutEntry(ptIDs, this->m_indexParams.r, mesh, rtree, this->m_indexParams.denseFlag, verbose); break;
		}
		return false;
	}

	void calRdBox(const float rd[RD_NUM], const float err[RD_NUM], int boxMin[INDEX_DIM], int boxMax[INDEX_DIM]);
	void calRdBox(const float rd[RD_NUM], int boxMin[INDEX_DIM], int boxMax[INDEX_DIM]);

	void analyzeRd(const float rd[RD_NUM], float& avgRd, float& minRd, float& maxRd, float& sdRd, float& sdRatioRd);

	const std::string getEntFilename(const std::string& path, const std::string& filename) {
		return fs::path(path).append(filename + ".ent").string();
	}

	const std::string getRstFilename(const std::string& path, const std::string& filename) {
		return fs::path(path).append(filename + ".rst").string();
	}


public:
	IndexParams m_indexParams;

	uint64_t m_numInserted;
	uint64_t m_numFailed;
	
	C2OIndex(MeshSet* _meshSet) {
		this->m_meshSet = _meshSet;
		this->m_numInserted = 0;
		this->m_numFailed = 0;
	}

	virtual ~C2OIndex() {
	}

	bool buildIndexForAll(const IndexParams& params, float* buildTime = nullptr, bool progressBar = false, bool verbose = false);
	
	bool buildIndexByPartsAndSave(const IndexParams& params, const std::string& path, const std::string& filename, float* stats = nullptr, bool progressBar = false, bool verbose = false);

	bool loadIndex(const std::string& path, const std::string& filename, bool verbose = false);

	bool saveIndex(const std::string& path, const std::string& filename, float* saveTime = nullptr, bool verbose = false);

	void analyzeIndexForAll(float& avg, float& min, float& max, float& sd, float& sdRatio);

	int searchEntries(const float rd[RD_NUM], const float err[RD_NUM], std::vector<const IndexEntry*>& retList, int* nodeAccessed = nullptr, int* leafAccessed = nullptr);
};

bool C2OIndex::buildIndexByPartsAndSave(const IndexParams& params, const std::string& path, const std::string& filename, float* stats, bool progressBar, bool verbose)
{
	if (this->c2oVariant != C2OV_DONUT) {
		std::cerr << "C2O variant not supported" << std::endl;
		return false;
	}
	if (params.r < 0.0f) {
		std::cerr << "Radius negative" << std::endl;
		return false;
	}

	this->m_indexParams = params;

	float avg = 0, min = 0, max = 0, sd = 0, sdRatio = 0;
	float avgRd, minRd, maxRd, sdRd, sdRatioRd;
	uint64_t analyzeRdCount = 0;

	float totalTime = 0, buildTime = 0, saveTime = 0, memoryTime = 0;

	if (stats) timer_start();

	if (stats) timer_start();

	const std::string entFilename = this->getEntFilename(path, filename);
	std::ofstream ofsEnt(entFilename, std::ios::out | std::ios::binary);
	if (!ofsEnt.is_open()) {
		std::cerr << "Error open " << entFilename << " for writing" << std::endl;
		return false;
	}

	const uint64_t total = this->m_meshSet->totalVertices();
	std::cout << "To Write " << total << " entries to " << entFilename << "..." << std::endl;
	
	ofsEnt.write((char *) &this->m_indexParams, sizeof(IndexParams));
	ofsEnt.write((char *) &total, sizeof(uint64_t));
	ofsEnt.flush();

	if (stats) saveTime += timer_end();


	uint64_t iGlobal = 0;
	int iPtIDs[4];
	bool entryCalculated;
	float rd[RD_NUM];
	int boxMin[INDEX_DIM], boxMax[INDEX_DIM];

	std::vector<IndexEntry*> entryListPart;
	int iPt;

    for (int iPart = 0; iPart < this->m_meshSet->getNumParts(); iPart++) {
    	if (stats) timer_start();
    	this->m_meshSet->buildRTreeForPart(iPart, progressBar, verbose);
    	if (stats) memoryTime += timer_end();

		if (verbose) {
			std::cout << "Building index for part #" << iPart << "... ";
			if (progressBar) std::cout << std::endl;
		}

    	entryListPart.clear();
    	std::vector<IndexEntry*>().swap(entryListPart);

    	ProgressBar bar(this->m_meshSet->getNumPtsForPart(iPart), 70);

    	for (int iMeshIdx = iPart * this->m_meshSet->getNumMeshesPerPart(); iMeshIdx < (iPart + 1) * this->m_meshSet->getNumMeshesPerPart(); iMeshIdx++) {
    		const auto& iMesh = this->m_meshSet->getMesh(iMeshIdx);

			for (iPt = 0; iPt < iMesh.mesh->vertices.size(); iPt++) {

		        if (progressBar) {
		            ++bar;
		            bar.display();
		        }

				iPtIDs[0] = iPt;

				// std::cout << "Global ID (" << iGlobal << "), Point ID (" << iPtIDs[0] << ")" << std::endl;

				entryCalculated = this->calEntry(iPtIDs, iMesh.mesh, iMesh.rtree);
				if (!entryCalculated) {
					for (int i = 1; i < 4; i++) {
						iPtIDs[i] = -1;
					}
				}

				// std::cout << "Returned IDs: " << iPtIDs[0] << ", " << iPtIDs[1] << ", " << iPtIDs[2] << ", " << iPtIDs[3] << std::endl;

				entryListPart.push_back(new IndexEntry(iMesh.id, iPtIDs));

				if (entryCalculated) { // calculate rd(side lengths and others) and insert into R*-tree
					calRd(iPtIDs, iMesh.mesh, rd);

					// for (int i = 0; i < RD_NUM; i++) {
					// 	std::cout << rd[i] << " ";
					// }
					// std::cout << std::endl;

					calRdBox(rd, boxMin, boxMax);

					// for (int i = 0; i < INDEX_DIM; i++) {
					// 	std::cout << "[" << boxMin[i] << ", " << boxMax[i] << "] ";
					// }
					// std::cout << std::endl;

					this->m_indexTree.Insert(boxMin, boxMax, iGlobal);
					
					this->m_numInserted++;

					#ifdef _IDX3
						for (int i = 0; i < 3; i++) {
							entryListPart.back()->unindexedSides[i] = rd[i];
						}

						// for (int i = 0; i < 3; i++) {
						// 	std::cout << this->m_entryList.back().unindexedSides[i] << " ";
						// }
						// std::cout << std::endl;
					#endif

					this->analyzeRd(rd, avgRd, minRd, maxRd, sdRd, sdRatioRd);

					avg += avgRd;
					min += minRd;
					max += maxRd;
					sd += sdRd;
					sdRatio += sdRatioRd;
					analyzeRdCount ++;

				} else {
					this->m_numFailed++;
				}

				iGlobal++;
			}
    	}

	    if (progressBar) {
	        bar.done();
	    }

		if (verbose) std::cout << "Done" << std::endl;

	    std::cout << "Writing part #" << iPart << " of " << entryListPart.size() << " entries to " << entFilename << "... ";
		if (stats) timer_start();
		for (uint64_t iEntry = 0; iEntry < entryListPart.size(); iEntry++) {
			ofsEnt.write((char *) (entryListPart[iEntry]), sizeof(IndexEntry)); // TODO: write a block
		}
		ofsEnt.flush();
		if (stats) saveTime += timer_end();
		std::cout << "Done" << std::endl;

    	if (stats) timer_start();
		for (uint64_t iEntry = 0; iEntry < entryListPart.size(); iEntry++) {
			delete entryListPart[iEntry];
		}
    	this->m_meshSet->freeRTreeForPart(iPart, verbose);
    	if (stats) memoryTime += timer_end();
    }

	ofsEnt.close();

    this->m_indexTree.SortDim0();

	const std::string rstFilename = this->getRstFilename(path, filename);

	std::cout << "Writing tree to " << rstFilename << "... ";
	if (stats) timer_start();
	this->m_indexTree.Save(rstFilename.c_str());
	if (stats) saveTime += timer_end();
	std::cout << "Done" << std::endl;

    if (stats) totalTime = timer_end();

    if (stats) {
    	buildTime = totalTime - memoryTime - saveTime;
    	stats[0] = buildTime;
    	stats[1] = saveTime;

		float fGlobalCount = float(analyzeRdCount);
		stats[2] = avg / fGlobalCount;
		stats[3] = min / fGlobalCount;
		stats[4] = max / fGlobalCount;
		stats[5] = sd / fGlobalCount;
		stats[6] = sdRatio / fGlobalCount;
    }

	if (verbose) {
		std::cout << "Total: " << total << std::endl;
		std::cout << "Inserted: " << this->m_numInserted << std::endl;
		std::cout << "Failed: " << this->m_numFailed << std::endl;
	}

	return true;
}

bool C2OIndex::buildIndexForAll(const IndexParams& params, float* buildTime, bool progressBar, bool verbose)
{
	if (this->c2oVariant != C2OV_DONUT) {
		std::cerr << "C2O variant not supported" << std::endl;
		return false;
	}
	if (params.r < 0.0f) {
		std::cerr << "Radius negative" << std::endl;
		return false;
	}

	this->m_indexParams = params;

	if (buildTime) timer_start();

	uint64_t iGlobal = 0;
	int iPtIDs[4];
	bool entryCalculated;
	float rd[RD_NUM];
	int boxMin[INDEX_DIM], boxMax[INDEX_DIM];

	const auto total = this->m_meshSet->totalVertices();
    
    ProgressBar bar(total, 70);

	for (const auto &iMesh: this->m_meshSet->getMeshList()) {
		for (int iPt = 0; iPt < iMesh.mesh->vertices.size(); iPt++) {

	        if (progressBar) {
	            ++bar;
	            bar.display();
	        }

			iPtIDs[0] = iPt;

			// std::cout << "Global ID (" << iGlobal << "), Point ID (" << iPtIDs[0] << ")" << std::endl;

			entryCalculated = this->calEntry(iPtIDs, iMesh.mesh, iMesh.rtree);
			if (!entryCalculated) {
				for (int i = 1; i < 4; i++) {
					iPtIDs[i] = -1;
				}
			}

			// std::cout << "Returned IDs: " << iPtIDs[0] << ", " << iPtIDs[1] << ", " << iPtIDs[2] << ", " << iPtIDs[3] << std::endl;

			this->m_entryList.emplace_back(iMesh.id, iPtIDs);

			if (entryCalculated) { // calculate rd(side lengths and others) and insert into R*-tree
				calRd(iPtIDs, iMesh.mesh, rd);

				// for (int i = 0; i < RD_NUM; i++) {
				// 	std::cout << rd[i] << " ";
				// }
				// std::cout << std::endl;

				calRdBox(rd, boxMin, boxMax);

				// for (int i = 0; i < INDEX_DIM; i++) {
				// 	std::cout << "[" << boxMin[i] << ", " << boxMax[i] << "] ";
				// }
				// std::cout << std::endl;

				this->m_indexTree.Insert(boxMin, boxMax, iGlobal);
				
				this->m_numInserted++;

				#ifdef _IDX3
					for (int i = 0; i < 3; i++) {
						this->m_entryList.back().unindexedSides[i] = rd[i];
					}

					// for (int i = 0; i < 3; i++) {
					// 	std::cout << this->m_entryList.back().unindexedSides[i] << " ";
					// }
					// std::cout << std::endl;
				#endif

			} else {
				this->m_numFailed++;
			}

			iGlobal++;
		}
	}

    if (progressBar) {
        bar.done();
    }

    this->m_indexTree.SortDim0();

    if (buildTime) (*buildTime) = timer_end();

	if (verbose) {
		std::cout << "Total: " << total << std::endl;
		std::cout << "Inserted: " << this->m_numInserted << std::endl;
		std::cout << "Failed: " << this->m_numFailed << std::endl;
	}

	return true;
}


bool C2OIndex::saveIndex(const std::string& path, const std::string& filename, float* saveTime, bool verbose)
{
	if (saveTime) timer_start();

	const std::string entFilename = this->getEntFilename(path, filename);
	std::ofstream ofsEnt(entFilename, std::ios::out | std::ios::binary);
	if (!ofsEnt.is_open()) {
		std::cerr << "Error open " << entFilename << " for writing" << std::endl;
		return false;
	}

	const uint64_t numEntries = this->m_entryList.size();
	std::cout << "Writing " << numEntries << " entries to " << entFilename << "... ";
	
	ofsEnt.write((char *) &this->m_indexParams, sizeof(IndexParams));
	ofsEnt.write((char *) &numEntries, sizeof(uint64_t));
	for (uint64_t iEntry = 0; iEntry < numEntries; iEntry++) {
		ofsEnt.write((char *) (&(this->m_entryList[iEntry])), sizeof(IndexEntry));
	}

	ofsEnt.close();
	std::cout << "Done" << std::endl;


	const std::string rstFilename = this->getRstFilename(path, filename);

	std::cout << "Writing tree to " << rstFilename << "... ";
	this->m_indexTree.Save(rstFilename.c_str());
	std::cout << "Done" << std::endl;

    if (saveTime) (*saveTime) = timer_end();

	return true;
}

bool C2OIndex::loadIndex(const std::string& path, const std::string& filename, bool verbose)
{
	const std::string entFilename = this->getEntFilename(path, filename);
	std::ifstream ifsEnt(entFilename, std::ios::in | std::ios::binary);
	if (!ifsEnt.is_open()) {
		std::cerr << "Error open " << entFilename << std::endl;
		return false;
	}

	if (verbose) {
		std::cout << "Reading entries from " << entFilename << "... ";
	}

	ifsEnt.read((char *) &this->m_indexParams, sizeof(IndexParams));

	uint64_t numEntries;
	ifsEnt.read((char *) &numEntries, sizeof(uint64_t));

	IndexEntry entry;
	for (uint64_t iEntry = 0; iEntry < numEntries; iEntry++) {
		ifsEnt.read((char *) &entry, sizeof(IndexEntry));
		this->m_entryList.push_back(entry);
	}

	ifsEnt.close();
	if (verbose) {
		std::cout << "Done (" << numEntries << ")" << std::endl;
	}

	const std::string rstFilename = this->getRstFilename(path, filename);

	if (verbose) {
		std::cout << "Reading tree from " << rstFilename << "... ";
	}
	this->m_indexTree.Load(rstFilename.c_str());
	if (verbose) {
		std::cout << "Done" << std::endl;
	}

	return true;
}


void C2OIndex::analyzeRd(const float rd[RD_NUM], float& avgRd, float& minRd, float& maxRd, float& sdRd, float& sdRatioRd)
{
	minRd = std::numeric_limits<float>::max();
	maxRd = -1;
	avgRd = 0.0;
	sdRd = 0.0;
	for (int i = 0; i < RD_NUM; i++) {
		avgRd += rd[i];
		if (maxRd < rd[i]) {
			maxRd = rd[i];
		}
		if (minRd > rd[i]) {
			minRd = rd[i];
		}
	}
	avgRd /= float(RD_NUM);
	for (int i = 0; i < RD_NUM; i++) {
		sdRd += (rd[i] - avgRd) * (rd[i] - avgRd);
	}
	sdRd /= float(RD_NUM);
	sdRd = sqrt(sdRd);
	sdRatioRd = sdRd / avgRd;
}


void C2OIndex::analyzeIndexForAll(float& avg, float& min, float& max, float& sd, float& sdRatio)
{
	avg = min = max = sd = sdRatio = 0;
	float rd[RD_NUM];
	uint64_t globalCount = 0;
	float avgRd, minRd, maxRd, sdRd, sdRatioRd;

	for (auto &entry: this->m_entryList) {
		if (!entry.isValid()) {
			continue;
		}

		calRd(entry.ptIDs, this->m_meshSet->getMesh(entry.dbID).mesh, rd);
		this->analyzeRd(rd, avgRd, minRd, maxRd, sdRd, sdRatioRd);

		avg += avgRd;
		min += minRd;
		max += maxRd;
		sd += sdRd;
		sdRatio += sdRatioRd;
		globalCount ++;
	}

	float fGlobalCount = float(globalCount);
	avg /= fGlobalCount;
	min /= fGlobalCount;
	max /= fGlobalCount;
	sd /= fGlobalCount;
	sdRatio /= fGlobalCount;
}

int C2OIndex::searchEntries(const float rd[RD_NUM], const float err[RD_NUM], std::vector<const IndexEntry*>& retList, int* nodeAccessed, int* leafAccessed)
{
	int boxMin[INDEX_DIM], boxMax[INDEX_DIM];
	calRdBox(rd, err, boxMin, boxMax);

	std::vector<int> treeRet;
	this->m_indexTree.Search(boxMin, boxMax, treeRet, nodeAccessed, leafAccessed);

	int num = 0;
	for (auto &iTreeRet: treeRet) {
		auto entry = (const IndexEntry*) (&(this->m_entryList[iTreeRet]));

		#ifdef _IDX3
			bool filtered = false;
			for (int i = 0; i < 3; i++) {
				if (abs(rd[i] - entry->unindexedSides[i]) > err[i]) {
					filtered = true;
					break;
				}
			}
			if (filtered) {
				continue;
			}
		#endif

		retList.push_back(entry);
		num++;
	}

	return num;
}

void C2OIndex::calRdBox(const float rd[RD_NUM], const float err[RD_NUM], int boxMin[INDEX_DIM], int boxMax[INDEX_DIM])
{
	#ifndef _CLR
		#ifdef _IDX3
			// use the last 3 side length as index keys
			boxMin[0] = int((std::max(0.0f, rd[3] - err[3])) * RSTREE_SCALE);
			boxMax[0] = int((rd[3] + err[3]) * RSTREE_SCALE);
			boxMin[1] = int((std::max(0.0f, rd[4] - err[4])) * RSTREE_SCALE);
			boxMax[1] = int((rd[4] + err[4]) * RSTREE_SCALE);
			boxMin[2] = int((std::max(0.0f, rd[5] - err[5])) * RSTREE_SCALE);
			boxMax[2] = int((rd[5] + err[5]) * RSTREE_SCALE);
		#else
			// use 6-side length as index keys
			for (int i = 0; i < INDEX_DIM; i++) {
				boxMin[i] = int((std::max(0.0f, rd[i] - err[i])) * RSTREE_SCALE);
				boxMax[i] = int((rd[i] + err[i]) * RSTREE_SCALE);
			}
		#endif
	#else
		// TODO
	#endif
}

void C2OIndex::calRdBox(const float rd[RD_NUM], int boxMin[INDEX_DIM], int boxMax[INDEX_DIM])
{
    float zeroErr[RD_NUM] = {};
    calRdBox(rd, zeroErr, boxMin, boxMax);
}


#endif