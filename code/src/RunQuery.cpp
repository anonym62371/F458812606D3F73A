#include <iostream>
#include <string>

#include "C2OIndex.hpp"
#include "C2OQuery.hpp"

using namespace std;
namespace fs = std::filesystem;
using namespace fs;

bool readConfigQuery(const string& configFile, QueryParams& params, string& queryPath, string& queryFilename, string& indexPath, string& indexFilename)
{
	ifstream ifs(configFile);
	if (!ifs.is_open()) {
		cout << "Error opening config file" << endl;
		return false;
	}

	string line;
	while (getline(ifs, line)) {
		if (isspace(line[line.length() - 1])) {
			line = line.substr(0, line.length() - 1);
		}
		auto line_cstr = line.c_str();

		if (line.rfind("delta=", 0) == 0) {
			params.delta = atof(line_cstr + 6);

		} else if (line.rfind("h=", 0) == 0) {
			params.hCF = atoi(line_cstr + 2);

		} else if (line.rfind("index_path=", 0) == 0) {
			indexPath = string(line_cstr + 11);

		} else if (line.rfind("index_filename=", 0) == 0) {
			indexFilename = string(line_cstr + 15);

		} else if (line.rfind("query_path=", 0) == 0) {
			queryPath = string(line_cstr + 11);

		} else if (line.rfind("query_filename=", 0) == 0) {
			queryFilename = string(line_cstr + 15);

		}
	}

	ifs.close();

	return true;
}

int main(int argc, char *argv[])
{
	trimesh::TriMesh::set_verbose(0);

	auto msDB = new MeshSet(string(argv[1]));
	msDB->buildKDTreeForAll();

	string configFilename = "config_query.txt";
	if (argc > 2) {
		configFilename = string(argv[2]);
	}

	QueryParams params;
	string queryPath = "", queryFilename = "";
	string indexPath = "", indexFilename = "";

	const string configFilePath = path(string(argv[1])).append(configFilename).string();
	if (!readConfigQuery(configFilePath, params, queryPath, queryFilename, indexPath, indexFilename)) {
		exit(1);
	}

	if (queryPath == "" || queryPath == "<default>") {
		queryPath = msDB->getRootPath();
	}
	if (queryFilename == "" || queryFilename == "<default>") {
		queryFilename = "query.ply";
	}

	if (indexPath == "" || indexPath == "<default>") {
		indexPath = msDB->getRootPath();
	}
	if (indexFilename == "" || indexFilename == "<default>") {
		indexFilename = "index";
	}

	const string queryFilePath = path(queryPath).append(queryFilename).string();
	cout << "Query file: " << queryFilePath << endl;
	cout << "Index path: " << indexPath << endl;
	cout << "Index filename: " << indexFilename << endl;
	cout << "delta: " << params.delta << endl;
	cout << "h: " << params.hCF << endl;

	auto c2oIndex = new C2OIndex(msDB);
	c2oIndex->loadIndex(indexPath, indexFilename);

	auto msQuery = new MeshSet(queryFilePath);
	msQuery->buildRTreeForAll();

	auto trimeshQuery = msQuery->getMesh(0).mesh;
	trimeshQuery->need_bsphere();
	float diamQuery = trimeshQuery->bsphere.r * 2.0;
	params.delta *= diamQuery;

	C2OQuery query(msDB, msQuery, c2oIndex);

    rusage start;
    getrusage(RUSAGE_SELF, &start);

    RetMatches rm;
	int matchedSize = query.execute(params, rm, nullptr, false);

	rusage end;
    getrusage(RUSAGE_SELF, &end);
    float queryTime = calculateExecutionTime(&start, &end);

	ofstream ofs(path(queryPath).append("output.txt").string());
	ofs << queryTime << endl;
	ofs << matchedSize << endl;

	for (int i = 0; i < matchedSize; i++) {
        ofs << rm[i].id << " " << msDB->getMesh(rm[i].id).filename << " " << rm[i].dist << endl;
        ofs << rm[i].xf;
	}

	cout << endl;

	delete msDB;
	delete c2oIndex;
	delete msQuery;
}