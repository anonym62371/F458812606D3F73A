#ifndef __STATSUTIL_H_
#define __STATSUTIL_H_

#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

class StatManager {
	std::vector<std::string> name_list;
	std::unordered_map<std::string, int> name_map;
	std::vector<std::string> value_list;
public:
	StatManager(const std::vector<std::string>& p_name_list) {
		name_list = p_name_list;
		for (int i = 0; i < name_list.size(); i++) {
			name_map[name_list[i]] = i;
			value_list.push_back("");
		}
	}
	StatManager() { }
	std::string get(const std::string& name) {
		if (name_map.find(name) == name_map.end()) {
			return "";
		} else {
			return value_list[name_map[name]];
		}
	}
	template <typename T>
	void set(const std::string& name, const T& t) {
		if (name_map.find(name) == name_map.end()) {
			name_list.push_back(name);
			name_map[name] = name_list.size() - 1;
			value_list.push_back("");
		}

		std::ostringstream os;
		os << t;
		value_list[name_map[name]] = os.str();
	}
	void write(std::ostream& os) {
		for (int i = 0; i < name_list.size(); i++) {
			os << name_list[i] << "=" << value_list[i] << std::endl;
		}
	}
};

#endif