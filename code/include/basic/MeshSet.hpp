#ifndef __MESHSET_H_
#define __MESHSET_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <filesystem>
#include <cstdint>
#include <cassert>
#include <cstdio>

#include "ProgressBar.hpp"
#include "C_RTree.hpp"
#include "UtilVec.hpp"
#include "UtilPolygon.hpp"
#include "UtilTrans.hpp"

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "KDtree.h"

namespace fs = std::filesystem;

struct SingleMesh
{
	uint32_t id; // the index of a mesh, might be as large as the total number of meshes
	std::string filename;
	trimesh::TriMesh* mesh;
	C_RTree* rtree;
	trimesh::KDtree* kdtree;
	std::vector<Triangle> triangles;
	std::vector<float> triangleAreas;
	float totalArea;
	std::vector<TriSide> edges;
	std::vector<std::pair<int, int>> edgeIDs;
	SingleMesh(uint32_t _id, const std::string& _filename, trimesh::TriMesh* _mesh) {
		this->id = _id;
		this->filename = _filename;
		this->mesh = _mesh;
		this->rtree = nullptr;
		this->kdtree = nullptr;
	}
};

class MeshSet
{
public:
	MeshSet() {
		this->numParts = 0;
		this->numMeshesPerPart = 0;

		this->rescaleRatio = 0;
	}

	MeshSet(const std::string& path, bool verbose = false, float rescaleRatio = 0) {
		this->numParts = 0;
		this->numMeshesPerPart = 0;

		this->rescaleRatio = rescaleRatio;

		read(path, verbose);
	}

	virtual ~MeshSet() {
		for (auto & iMesh: this->meshList) {
			if (iMesh.mesh) {
				iMesh.mesh->clear();
				delete iMesh.mesh;
				iMesh.mesh = nullptr;
			}
			if (iMesh.rtree) {
				delete iMesh.rtree;
				iMesh.rtree = nullptr;
			}
			if (iMesh.kdtree) {
				delete iMesh.kdtree;
				iMesh.kdtree = nullptr;
			}
		}
		std::vector<SingleMesh>().swap(meshList);
	}

	uint32_t read(const std::string& path, bool verbose = false);

	void buildRTreeForAll(bool verbose = false);
	void saveRTreeForAll(bool verbose = false);

	bool buildRTreeForPart(int iPart, bool progressBar = false, bool verbose = false);
	bool freeRTreeForPart(int iPart, bool verbose = false);

	void buildKDTreeForAll(bool verbose = false);
	
	float calDist2(const trimesh::TriMesh* mesh, uint32_t dbID, const trimesh::xform* xf = nullptr, float stopAtErr2 = std::numeric_limits<float>::max());

	void loadTrianglesForAll(bool verbose = false);
	void loadTriangleAreasForAll(bool verbose = false);

	void loadEdgesForAll(bool verbose = false);

	void loadBBoxForAll(bool verbose = false);
	void loadBSphereForAll(bool verbose = false);

	void loadNormalForAll(bool verbose = false);

	const std::vector<SingleMesh> & getMeshList() {
		return this->meshList;
	}

	SingleMesh & getMesh(int _index) {
		assert(_index >= 0 && _index < this->meshList.size());
		return this->meshList[_index];
	}

	void addMesh(const SingleMesh & _mesh) {
		this->meshList.push_back(_mesh);
		this->totalV += _mesh.mesh->vertices.size();
	}

	uint32_t size() const {
		return this->meshList.size();
	}
	uint64_t totalVertices() const {
		return this->totalV;
	}
	uint64_t totalFaces() const {
		return this->totalF;
	}

	std::string getRootPath() const {
		return rootPath.string();
	}

	uint32_t getNumParts() const {
		return numParts;
	}

	uint32_t getNumMeshesPerPart() const {
		return numMeshesPerPart;
	}

	uint64_t getNumPtsForPart(int iPart) const {
		if (iPart < 0 || iPart >= this->numParts) {
			return 0;
		}
		return numPtsEachPart[iPart];
	}

	void setNumParts(uint32_t _numParts) {
		this->numParts = _numParts;
	}

	// void save() const; // TODO

private:
	std::vector<SingleMesh> meshList;
	uint64_t totalV;
	uint64_t totalF;

	fs::path rootPath;
	uint32_t numParts;
	uint32_t numMeshesPerPart;
	std::vector<uint64_t> numPtsEachPart;

	float rescaleRatio;

	trimesh::TriMesh* readTriMesh(const std::string & filename) const;

	bool readMeta(uint32_t& numMeshes, std::vector<uint32_t>& meshIDList, std::vector<std::string>& meshNameList, bool verbose);
	uint32_t readDir(bool verbose);
	uint32_t readComb(bool verbose);

	trimesh::TriMesh* meshSelect(const trimesh::TriMesh* meshFrom, int64_t ptFrom = -1, int64_t ptSize = -1) const;
	const std::string getPath(const std::string& filename) const;

	bool existsInEdgeHash(const std::unordered_map<int, std::unordered_set<int>*> & edgeHash, int ep1ID, int ep2ID);
	void insertIntoEdgeHash(std::unordered_map<int, std::unordered_set<int>*> & edgeHash, int ep1ID, int ep2ID);
	void deleteEdgeHash(std::unordered_map<int, std::unordered_set<int>*> & edgeHash);

};


void MeshSet::loadTrianglesForAll(bool verbose)
{
	if (verbose) {
		std::cout << "Loading triangles for all... ";
	}

	this->totalF = 0;

	for (auto & iMesh: this->meshList) {
		iMesh.mesh->need_faces();
		iMesh.mesh->need_adjacentfaces();
		for (int iFace = 0; iFace < iMesh.mesh->faces.size(); iFace++) {
			iMesh.triangles.emplace_back(iFace, iMesh.mesh);

			this->totalF ++;
		}
	}

	if (verbose) {
		std::cout << "Done" << std::endl;
	}
}

void MeshSet::loadTriangleAreasForAll(bool verbose)
{
	if (verbose) {
		std::cout << "Loading triangle areas for all... ";
	}

	for (auto & iMesh: this->meshList) {
		iMesh.totalArea = 0;
		for (auto & iTri: iMesh.triangles) {
			float area = areaTri(iTri);
			iMesh.triangleAreas.push_back(area);
			iMesh.totalArea += area;
		}
	}

	if (verbose) {
		std::cout << "Done" << std::endl;
	}
}

bool MeshSet::existsInEdgeHash(const std::unordered_map<int, std::unordered_set<int>*> & edgeHash, int ep1ID, int ep2ID)
{
	int minID = std::min(ep1ID, ep2ID);
	int maxID = std::max(ep1ID, ep2ID);

	if (edgeHash.find(minID) == edgeHash.end()) {
		return false;
	}

	if (edgeHash.at(minID)->find(maxID) == edgeHash.at(minID)->end()) {
		return false;
	} else {
		return true;
	}
}

void MeshSet::insertIntoEdgeHash(std::unordered_map<int, std::unordered_set<int>*> & edgeHash, int ep1ID, int ep2ID)
{
	int minID = std::min(ep1ID, ep2ID);
	int maxID = std::max(ep1ID, ep2ID);

	if (edgeHash.find(minID) == edgeHash.end()) {
		edgeHash[minID] = new std::unordered_set<int>;
	}
	edgeHash[minID]->insert(maxID);
}

void MeshSet::deleteEdgeHash(std::unordered_map<int, std::unordered_set<int>*> & edgeHash)
{
	for (auto & p: edgeHash) {
		delete p.second;
	}
	edgeHash.clear();
}

void MeshSet::loadEdgesForAll(bool verbose)
{
	if (verbose) {
		std::cout << "Loading edges for all... ";
	}

	for (auto & iMesh: this->meshList) {
		iMesh.mesh->need_faces();
		std::unordered_map<int, std::unordered_set<int>*> edgeHash; // TODO: change into number hash?
		
		for (auto & iFace: iMesh.mesh->faces) {
			if (!this->existsInEdgeHash(edgeHash, iFace[0], iFace[1])) {
				iMesh.edges.emplace_back(iFace[0], iFace[1], iMesh.mesh);
				iMesh.edgeIDs.emplace_back(iFace[0], iFace[1]);
				this->insertIntoEdgeHash(edgeHash, iFace[0], iFace[1]);
			}
			if (!this->existsInEdgeHash(edgeHash, iFace[0], iFace[2])) {
				iMesh.edges.emplace_back(iFace[0], iFace[2], iMesh.mesh);
				iMesh.edgeIDs.emplace_back(iFace[0], iFace[2]);
				this->insertIntoEdgeHash(edgeHash, iFace[0], iFace[2]);
			}
			if (!this->existsInEdgeHash(edgeHash, iFace[1], iFace[2])) {
				iMesh.edges.emplace_back(iFace[1], iFace[2], iMesh.mesh);
				iMesh.edgeIDs.emplace_back(iFace[1], iFace[2]);
				this->insertIntoEdgeHash(edgeHash, iFace[1], iFace[2]);
			}
		}

		this->deleteEdgeHash(edgeHash);
	}

	if (verbose) {
		std::cout << "Done" << std::endl;
	}
}

void MeshSet::loadBBoxForAll(bool verbose)
{
	if (verbose) {
		std::cout << "Loading bounding boxes for all... ";
	}

	for (auto & iMesh: this->meshList) {
		iMesh.mesh->need_bbox();
	}

	if (verbose) {
		std::cout << "Done" << std::endl;
	}
}

void MeshSet::loadBSphereForAll(bool verbose)
{
	if (verbose) {
		std::cout << "Loading bounding spheres for all... ";
	}

	for (auto & iMesh: this->meshList) {
		iMesh.mesh->need_bsphere();
	}

	if (verbose) {
		std::cout << "Done" << std::endl;
	}
}

void MeshSet::loadNormalForAll(bool verbose)
{
	if (verbose) {
		std::cout << "Loading normals for all... ";
	}

	for (auto & iMesh: this->meshList) {
		iMesh.mesh->need_normals();
	}

	if (verbose) {
		std::cout << "Done" << std::endl;
	}
}

void MeshSet::buildKDTreeForAll(bool verbose)
{
	if (verbose) {
		std::cout << "Building KD-tree for all... ";
	}

	for (auto & iMesh: this->meshList) {
		iMesh.kdtree = new trimesh::KDtree(iMesh.mesh->vertices);
	}

	if (verbose) {
		std::cout << "Done" << std::endl;
	}
}

float MeshSet::calDist2(const trimesh::TriMesh* mesh, uint32_t dbID, const trimesh::xform* xf, float stopAtErr2)
{
	if (!meshList[dbID].kdtree) {
		std::cerr << "KD-tree not built for mesh #" << dbID << std::endl;
		return -1.0f;
	}
    float err2 = 0.0f;
    trimesh::point xf_v;
    for (auto &v: mesh->vertices) {
    	if (xf) {
        	trans_pt(*xf, v, xf_v);
        } else {
        	xf_v = v;
        }

        // Original code: using kd-tree
        // auto nn = meshList[dbID].kdtree->closest_to_pt(xf_v);
        // add this big float to avoid the issue that some points may not find its NN because it is too far away
        // TODO: not verifying its performance 
        auto nn = meshList[dbID].kdtree->closest_to_pt(xf_v, 3.3e33);
        if (!nn) {
        	std::cout << "a point has null-NN for dbID: " << dbID << std::endl;
            continue;
        }
        auto nnPt = (const trimesh::point *) nn;
        float nnDist2 = dist2(*nnPt, xf_v);

        // // tmp code: bf-finding the nn
        // float nnDist2 = 3.3e33;
        // for (auto &mesh_v: meshList[dbID].mesh->vertices) {
        // 	float thisDist2 = dist2(mesh_v, xf_v);
        // 	if (nnDist2 > thisDist2) {
        // 		nnDist2 = thisDist2;
        // 	}
        // }
        
        err2 += nnDist2;
        if (err2 > stopAtErr2) {
            return -1.0f;
        }
    }
    return err2;
}

bool MeshSet::buildRTreeForPart(int iPart, bool progressBar, bool verbose)
{
	if (iPart < 0 || iPart >= this->numParts) {
		return false;
	}

	auto rstPath = this->getPath("combined.ply." + std::to_string(iPart) + ".rst.0");
	FILE *rstFilePtr = fopen(rstPath.c_str(), "r");
	if (rstFilePtr) { // if exists the r-tree file for the part
		if (verbose) {
			std::cout << "Loading R-tree for part " << iPart << " from " << rstPath << "... " << std::flush;
		}

		for (uint32_t iMesh = iPart * this->numMeshesPerPart; iMesh < (iPart + 1) * this->numMeshesPerPart; iMesh++) {
			this->meshList[iMesh].rtree = new C_RTree;
			this->meshList[iMesh].rtree->read(rstFilePtr);
			// std::cout << "R-tree loaded for mesh " << iMesh << std::endl;
		}

		if (verbose) {
			std::cout << "Done" << std::endl;
		}

		fclose(rstFilePtr);
	} else { // otherwise rebuild this part
		if (verbose) {
			std::cout << "Building R-tree for part " << iPart << "... ";
			if (progressBar) {
				std::cout << std::endl;
			}
		}

    	ProgressBar bar(this->numMeshesPerPart, 70);

		for (uint32_t iMesh = iPart * this->numMeshesPerPart; iMesh < (iPart + 1) * this->numMeshesPerPart; iMesh++) {
	        if (progressBar) {
	            ++bar;
	            bar.display();
	        }

			this->meshList[iMesh].rtree = new C_RTree;
			this->meshList[iMesh].rtree->build(this->meshList[iMesh].mesh);
			// std::cout << "R-tree built for mesh " << iMesh << std::endl;
		}

	    if (progressBar) {
	        bar.done();
	    }
		
		if (verbose) {
			std::cout << "Done" << std::endl;
		}
	}

	return true;
}

bool MeshSet::freeRTreeForPart(int iPart, bool verbose)
{
	if (iPart < 0 || iPart >= this->numParts) {
		return false;
	}

	if (verbose) {
		std::cout << "Free R-tree for part " << iPart << "... " << std::flush;
	}

	for (uint32_t iMesh = iPart * this->numMeshesPerPart; iMesh < (iPart + 1) * this->numMeshesPerPart; iMesh++) {
		delete this->meshList[iMesh].rtree;
		this->meshList[iMesh].rtree = nullptr;
		// std::cout << "R-tree deleted for mesh " << iMesh << std::endl;
	}
		
	if (verbose) {
		std::cout << "Done" << std::endl;
	}

	return true;
}

void MeshSet::buildRTreeForAll(bool verbose)
{
	if (this->numParts > 0) {
		// try to find saved r-tree by part
		const uint32_t numMeshes = this->meshList.size();
		uint32_t iGlobal = 0;
		for (uint32_t iPart = 0; iPart < this->numParts; iPart++) {
			if (!(iGlobal < numMeshes)) {
				break;
			}

			auto rstPath = this->getPath("combined.ply." + std::to_string(iPart) + ".rst.0");
			FILE *rstFilePtr = fopen(rstPath.c_str(), "r");
			if (rstFilePtr) { // if exists the r-tree file for the part
				if (verbose) {
					std::cout << "Loading R-tree for part " << iPart << " from " << rstPath << "... " << std::flush;
				}

				for (uint32_t iPly = 0; iPly < this->numMeshesPerPart; iPly++) {
					if (!(iGlobal < numMeshes)) {
						break;
					}
					this->meshList[iGlobal].rtree = new C_RTree;
					this->meshList[iGlobal].rtree->read(rstFilePtr);

					// std::cout << "R-tree loaded for mesh " << iGlobal << std::endl;

					iGlobal++;
				}

				if (verbose) {
					std::cout << "Done" << std::endl;
				}

				fclose(rstFilePtr);
			} else { // otherwise rebuild this part
				if (verbose) {
					std::cout << "Building R-tree for part " << iPart << "... " << std::flush;
				}

				for (uint32_t iPly = 0; iPly < this->numMeshesPerPart; iPly++) {
					if (!(iGlobal < numMeshes)) {
						break;
					}
					this->meshList[iGlobal].rtree = new C_RTree;
					this->meshList[iGlobal].rtree->build(this->meshList[iGlobal].mesh);

					// std::cout << "R-tree built for mesh " << iGlobal << std::endl;

					iGlobal++;
				}
				
				if (verbose) {
					std::cout << "Done" << std::endl;
				}
			}
		}
	} else {
		// try to find saved r-tree individually
		for (auto & iMesh: this->meshList) {
			iMesh.rtree = new C_RTree;
			auto rstPath = this->getPath(iMesh.filename + ".rst.0");
			FILE *rstFilePtr = fopen(rstPath.c_str(), "r");

			if (rstFilePtr) { // if exists the r-tree file ending up with ".rst.0"
				if (verbose) {
					std::cout << "Loading R-tree for " << iMesh.filename << " from " << rstPath << "... " << std::flush;
				}

				iMesh.rtree->read(rstFilePtr);
				
				if (verbose) {
					std::cout << "Done" << std::endl;
				}

				fclose(rstFilePtr);
			} else { // otherwise, rebuild it
				if (verbose) {
					std::cout << "Building R-tree for " << iMesh.filename << "... " << std::flush;
				}

				iMesh.rtree->build(iMesh.mesh);
				
				if (verbose) {
					std::cout << "Done" << std::endl;
				}
			}
		}
	}
}

// C++17 required
uint32_t MeshSet::read(const std::string& path, bool verbose)
{
	uint32_t readSize = 0;
	fs::path fsPath(path);

	if (fs::is_directory(fsPath)) {
		if (verbose) {
			std::cout << "Read directory " << fsPath << std::endl;
		}

		this->rootPath = fsPath;

		readSize = this->readDir(verbose);

	} else {
		if (verbose) {
			std::cout << "Read file " << fsPath << std::endl;
		}

		this->rootPath = fsPath.parent_path();

		// auto mesh = trimesh::TriMesh::read(fsPath.string());
		auto mesh = this->readTriMesh(fsPath.string());
		if (mesh) {
			this->meshList.emplace_back(0, fsPath.filename(), mesh);
			readSize = 1;
		}
	}

	if (verbose) {
		std::cout << "Read size " << readSize << std::endl;
		// std::cout << "Root path: " << rootPath << std::endl;
	}

	this->totalV = 0;
	for (const auto &iMesh: this->meshList) {
		this->totalV += iMesh.mesh->vertices.size();
	}
	if (verbose) {
		std::cout << "Total vertices: " << this->totalV << std::endl;
	}

	return readSize;
}

// C++17 required
bool MeshSet::readMeta(uint32_t& numMeshes, std::vector<uint32_t>& meshIDList, std::vector<std::string>& meshNameList, bool verbose)
{
	std::ifstream ifsMeta(this->getPath("meta.txt"));
	if (!ifsMeta.is_open()) {
		return false;
	}

	ifsMeta >> numMeshes;
	for (uint32_t iMesh = 0; iMesh < numMeshes; iMesh++) {
		uint32_t id;
		std::string name;
		ifsMeta >> id >> name;
		meshIDList.push_back(id);
		meshNameList.push_back(name);
	}

	ifsMeta.close();

	if (verbose) {
		std::cout << "Read meta of size " << numMeshes << std::endl;
	}

	return true;
}

// C++17 required
uint32_t MeshSet::readDir(bool verbose)
{
	// 1. If there exists combined.meta, then read combined
	std::ifstream ifsCombMeta(this->getPath("combined.meta"));
	if (ifsCombMeta.is_open()) {
		ifsCombMeta >> this->numParts >> this->numMeshesPerPart;
		ifsCombMeta.close();

		if (verbose) {
			std::cout << "Read combined in " << this->numParts << " parts each of size " << this->numMeshesPerPart << std::endl;
		}
		return this->readComb(verbose);
	}

	// 2. Otherwise, if meta.txt exists, read files according to meta
	uint32_t numMeshes = 0;
	uint32_t numMeshesReal = 0;
	std::vector<uint32_t> meshIDList;
	std::vector<std::string> meshNameList;

	if (this->readMeta(numMeshes, meshIDList, meshNameList, verbose)) {

		for (uint32_t iMesh = 0; iMesh < numMeshes; iMesh++) {

			// auto mesh = trimesh::TriMesh::read(this->getPath(meshNameList[iMesh]));
			auto mesh = this->readTriMesh(this->getPath(meshNameList[iMesh]));
			if (mesh) {
				this->meshList.emplace_back(meshIDList[iMesh], meshNameList[iMesh], mesh);
				numMeshesReal++;
			}

		}

		return numMeshesReal;
	}

	// 3. Otherwise, search for the files ending up with .ply
	uint32_t numPlyFound = 0;
    for (const auto &entry: fs::directory_iterator(rootPath)) {
        auto entryPath = entry.path();
        // std::cout << entryPath << std::endl;
        if (fs::is_directory(entryPath) || (entryPath.extension() != ".ply")) {
        	continue;
        }

        // auto mesh = trimesh::TriMesh::read(entryPath.string());
        auto mesh = this->readTriMesh(entryPath.string());
        if (mesh) {
        	this->meshList.emplace_back(numPlyFound, entryPath.filename(), mesh);
        	numPlyFound++;
        }

    }

    return numPlyFound;
}

// C++17 required
uint32_t MeshSet::readComb(bool verbose)
{
	uint32_t numMeshes;
	std::vector<uint32_t> meshIDList;
	std::vector<std::string> meshNameList;

	auto metaExists = this->readMeta(numMeshes, meshIDList, meshNameList, verbose);
	if (!metaExists) {
		numMeshes = this->numParts * this->numMeshesPerPart;
		for (uint32_t iMesh = 0; iMesh < numMeshes; iMesh++) {
			meshIDList.push_back(iMesh);
			meshNameList.push_back("mesh_" + std::to_string(iMesh));
		}
	}

	uint32_t iMesh = 0;
	bool incompleteRead = false;
	for (uint32_t iPart = 0; iPart < this->numParts; iPart++) {	
		if (!(iMesh < numMeshes)) {
			break;
		}

		std::vector<int64_t> plySizeList(this->numMeshesPerPart, 0);

		auto iPathCombMeta = this->getPath("combined.meta." + std::to_string(iPart));
		std::ifstream iIfsCombMeta(iPathCombMeta);
		if (iIfsCombMeta.is_open()) {
			if (verbose) {
				std::cout << "Read " << iPathCombMeta << std::endl;
			}

			for (uint32_t iPly = 0; iPly < this->numMeshesPerPart; iPly++) {
				iIfsCombMeta >> plySizeList[iPly];
			}

			iIfsCombMeta.close();
		} else {
			if (verbose) {
				std::cerr << "Read " << iPathCombMeta << " error" << std::endl;
			}
			incompleteRead = true;
			break;
		}

		// trimesh::TriMesh* iMeshComb = trimesh::TriMesh::read(this->getPath("combined.ply." + std::to_string(iPart)));
		trimesh::TriMesh* iMeshComb = this->readTriMesh(this->getPath("combined.ply." + std::to_string(iPart)));
		if (!iMeshComb) {
			incompleteRead = true;
			break;
		}

		int64_t ptFrom = 0;
		for (uint32_t iPly = 0; iPly < this->numMeshesPerPart; iPly++) {
			if (!(iMesh < numMeshes)) {
				break;
			}

			this->meshList.emplace_back(meshIDList[iMesh], meshNameList[iMesh], this->meshSelect(iMeshComb, ptFrom, plySizeList[iPly]));

			ptFrom += plySizeList[iPly];
			iMesh++;
		}

		this->numPtsEachPart.push_back(ptFrom);

        iMeshComb->clear();
        delete iMeshComb;

        std::vector<int64_t>().swap(plySizeList);
	}

	// std::cout << "Num pts each part:";
	// for (auto &v: this->numPtsEachPart) {
	// 	std::cout << " " << v;
	// }
	// std::cout << std::endl;

	if (incompleteRead) {
		return 0;
	} else {
		return numMeshes;
	}
}

// C++17 required
const std::string MeshSet::getPath(const std::string& filename) const
{
	fs::path fsPathConcat(rootPath);
	fsPathConcat /= filename;
	return fsPathConcat.string();
}

// TODO: future support: select faces
trimesh::TriMesh* MeshSet::meshSelect(const trimesh::TriMesh* meshFrom, int64_t ptFrom, int64_t ptSize) const
{
	if (ptFrom < 0) {
		ptFrom = 0;
	}
	if (ptSize < 0) {
		ptSize = meshFrom->vertices.size();
	}

	assert(ptFrom < meshFrom->vertices.size());
	assert((ptFrom + ptSize) <= meshFrom->vertices.size());

	trimesh::TriMesh* meshRet = new trimesh::TriMesh;
	for (int64_t i = 0; i < ptSize; i++) {
		meshRet->vertices.push_back(meshFrom->vertices[ptFrom + i]);
	}

	return meshRet;
}

trimesh::TriMesh* MeshSet::readTriMesh(const std::string & filename) const
{
	trimesh::TriMesh* mesh = trimesh::TriMesh::read(filename);
	if (this->rescaleRatio != 0) {
		scale(mesh, this->rescaleRatio);
	}
	// mesh->need_bsphere();
	// std::cout << "Mesh diam: " << (mesh->bsphere.r * 2) << std::endl;
	return mesh;
}

#endif