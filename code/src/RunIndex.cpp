#include <iostream>
#include <string>
#include <filesystem>

#include "C2OIndex.hpp"

using namespace std;
namespace fs = std::filesystem;
using namespace fs;

bool readConfigIndex(const string& configFile, IndexParams& params, string& indexPath, string& indexFilename, bool& showProgress)
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

		if (line.rfind("r=", 0) == 0) {
			params.r = atof(line_cstr + 2);

		} else if (line.rfind("index_path=", 0) == 0) {
			indexPath = string(line_cstr + 11);

		} else if (line.rfind("index_filename=", 0) == 0) {
			indexFilename = string(line_cstr + 15);

		} else if (line.rfind("show_progress_bar=", 0) == 0) {
			showProgress = atoi(line_cstr + 18);
		
		}
	}

	ifs.close();

	return true;
}

int main(int argc, char *argv[])
{
	auto ms = new MeshSet(string(argv[1]));
	ms->buildRTreeForAll();

	string configFilename = "config_index.txt";
	if (argc > 2) {
		configFilename = string(argv[2]);
	}

	IndexParams params;
	string indexPath = "", indexFilename = "";
	bool showProgress;

	const string configFilePath = path(string(argv[1])).append(configFilename).string();
	if (!readConfigIndex(configFilePath, params, indexPath, indexFilename, showProgress)) {
		exit(1);
	}

	if (indexPath == "" || indexPath == "<default>") {
		indexPath = ms->getRootPath();
	}
	if (indexFilename == "" || indexFilename == "<default>") {
		indexFilename = "index";
	}

	printIndexParams(params);
	cout << "Index path: " << indexPath << endl;
	cout << "Index filename: " << indexFilename << endl;
	// cout << "Show progress bar: " << showProgress << endl;

	C2OIndex c2oIndex(ms);
	c2oIndex.buildIndexForAll(params, nullptr, showProgress);
	c2oIndex.saveIndex(indexPath, indexFilename);

	cout << endl;

	C2OIndex c2oIndexLoad(ms);
	c2oIndexLoad.loadIndex(indexPath, indexFilename);

	cout << endl;

	delete ms;
}