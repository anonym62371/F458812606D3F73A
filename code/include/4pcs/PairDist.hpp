#ifndef __PAIRDIST_H_
#define __PAIRDIST_H_

#include <iostream>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <fstream>

#include "TriMesh.h"

#include "ProgressBar.hpp"
#include "UtilTime.hpp"
#include "MeshSet.hpp"

namespace fs = std::filesystem;

struct ParamsPD {
	float binSize;
};

void printParamsPD(const ParamsPD& params)
{
	std::cout << "Parameters:" << std::endl;
	std::cout << "bin-size: " << params.binSize << std::endl;
}

struct PairDistEntry
{
	int m_db_id;
	int m_ep1;
	int m_ep2;
	float m_dist;
	PairDistEntry() { }
	PairDistEntry(int db_id, int ep1, int ep2, float dist): m_db_id {db_id}, m_ep1 {ep1}, m_ep2 {ep2}, m_dist {dist} { }
};

struct PairDistBin
{
	int m_bin_id;
	std::vector<PairDistEntry> m_pd_list;
	PairDistBin() {}
	PairDistBin(int bin_id): m_bin_id {bin_id} {}
	void insert(int db_id, int ep1, int ep2, float dist) {
		m_pd_list.emplace_back(db_id, ep1, ep2, dist);
	}
	struct BinStruct {
		int bin_id;
		int pd_list_size;
	};
};

class PairDist
{
private:
	MeshSet* m_meshSet;

	std::unordered_map<int, PairDistBin> m_bin_map;
	long long m_total_pd;

	void insert(int db_id, int ep1, int ep2, float dist);

	int get_bin_id(float dist) {
		return static_cast<int>(dist / m_paramsPD.binSize);
	}

	struct BinMetaStruct {
		float binSize;
		int mapSize;
	};

public:
	ParamsPD m_paramsPD;
	
	PairDist(MeshSet* _meshSet) {
		m_meshSet = _meshSet;
		m_total_pd = 0;
	}

	virtual ~PairDist() {
	}

	int getNumBins() {
		return m_bin_map.size();
	}

	long long getTotalPD() {
		return m_total_pd;
	}

	bool buildIndexForAll(const ParamsPD& params, float* buildTime = nullptr, bool progressBar = false, bool verbose = false);

	bool saveIndex(const fs::path& path, const std::string& filename, float* saveTime = nullptr, bool verbose = false);

	bool loadIndex(const fs::path& path, const std::string& filename, bool verbose = false);

	void searchEntries(float dist_from, float dist_to, std::vector<PairDistEntry*>& ret);
};

void PairDist::insert(int db_id, int ep1, int ep2, float dist)
{
	int bin_id = this->get_bin_id(dist);
	if (m_bin_map.find(bin_id) == m_bin_map.end()) {
		m_bin_map[bin_id] = PairDistBin(bin_id);
	}
	m_bin_map[bin_id].insert(db_id, ep1, ep2, dist);
}

bool PairDist::buildIndexForAll(const ParamsPD& params, float* buildTime, bool progressBar, bool verbose)
{
	if (params.binSize < 0.0f) {
		std::cerr << "bin-size negative" << std::endl;
		return false;
	}

	this->m_paramsPD = params;

	if (buildTime) timer_start();

    trimesh::TriMesh* mesh;
    int meshSize;
    long long totalPairs;
    ProgressBar* bar;
    int iPt, jPt;

	for (int iMesh = 0; iMesh < this->m_meshSet->size(); iMesh++) {

		mesh = this->m_meshSet->getMesh(iMesh).mesh;
		meshSize = mesh->vertices.size();
		totalPairs = ((long long) meshSize) * (meshSize - 1) / 2;

        if (progressBar) {    
            bar = new ProgressBar(totalPairs, 70);
        }

		for (iPt = 0; iPt < meshSize; iPt++) {
			for (jPt = iPt + 1; jPt < meshSize; jPt++) {

		        if (progressBar) {
		            ++(*bar);
		            bar->display();
		        }

		        this->insert(iMesh, iPt, jPt, dist(mesh->vertices[iPt], mesh->vertices[jPt]));

		        this->m_total_pd ++;

		    }
		}

		if (progressBar) {
            bar->done();
            delete bar;
		}
	}

    if (buildTime) (*buildTime) = timer_end();

    return true;
}

bool PairDist::saveIndex(const fs::path& path, const std::string& filename, float* saveTime, bool verbose)
{
	if (saveTime) timer_start();

	const auto bin_filename = path / (filename + ".bin");
	std::ofstream ofs(bin_filename, std::ios::out | std::ios::binary);
	if (!ofs.is_open()) {
		std::cerr << "Fail open " << bin_filename << " for write" << std::endl;
	}

	// cout << sizeof(BinMetaStruct) << endl;
	// cout << sizeof(PairDistBin::BinStruct) << endl;
	// cout << sizeof(PairDist) << endl;

	BinMetaStruct bms { m_paramsPD.binSize, (int) m_bin_map.size() };
	ofs.write((char *) &bms, sizeof(BinMetaStruct));

	for (auto &v: m_bin_map) {
		auto bin = v.second;
		PairDistBin::BinStruct pdb_bs { bin.m_bin_id, (int) bin.m_pd_list.size() };
		ofs.write((char *) &pdb_bs, sizeof(PairDistBin::BinStruct));
		for (auto &p: bin.m_pd_list) {
			ofs.write((char *) &p, sizeof(PairDistEntry));
		}
	}

	ofs.close();

    if (saveTime) (*saveTime) = timer_end();

	return true;
}

bool PairDist::loadIndex(const fs::path& path, const std::string& filename, bool verbose)
{
	const auto bin_filename = path / (filename + ".bin");
	std::ifstream ifs(bin_filename, std::ios::in | std::ios::binary);
	if (!ifs.is_open()) {
		std::cerr << "Fail open " << bin_filename << " for read" << std::endl;
		return false;
	}

	BinMetaStruct bms;
	ifs.read((char *) &bms, sizeof(BinMetaStruct));
	this->m_paramsPD.binSize = bms.binSize;

	for (int i = 0; i < bms.mapSize; i++) {
		PairDistBin::BinStruct pdb_bs;
		ifs.read((char *) &pdb_bs, sizeof(PairDistBin::BinStruct));

		PairDistBin pdb(pdb_bs.bin_id);
		for (int j = 0; j < pdb_bs.pd_list_size; j++) {
			PairDistEntry pd;
			ifs.read((char *) &pd, sizeof(PairDistEntry));
			pdb.m_pd_list.push_back(pd);
		}

		m_bin_map[pdb_bs.bin_id] = pdb;
	}

	return true;
}

void PairDist::searchEntries(float dist_from, float dist_to, std::vector<PairDistEntry*>& ret)
{
	auto id_from = this->get_bin_id(dist_from);
	auto id_to = this->get_bin_id(dist_to);

	for (int id_i = id_from; id_i <= id_to; id_i++) {
		auto ptr_bin = &this->m_bin_map[id_i];
		if (id_i != id_from && id_i != id_to) {
			for (int i = 0; i < ptr_bin->m_pd_list.size(); i++) {
				ret.push_back(&ptr_bin->m_pd_list[i]);
			}
		} else {
			for (int i = 0; i < ptr_bin->m_pd_list.size(); i++) {
				if (dist_from <= ptr_bin->m_pd_list[i].m_dist && ptr_bin->m_pd_list[i].m_dist <= dist_to) {
					ret.push_back(&ptr_bin->m_pd_list[i]);
				}
			}
		}
	}
}


#endif