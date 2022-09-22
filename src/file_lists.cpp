#include <iostream>
#include <boost/filesystem.hpp>
#include <sstream>
#include <boost/regex.hpp>
#include "file_lists.hpp"
#include <algorithm>

using namespace std;
using namespace boost::filesystem;

string get_dat_name(string name, string dir){
	ostringstream full_name;
	full_name << dir << name;
	return full_name.str();
}

string get_log_name(string base, string dir, string tag){
        ostringstream log_stream;

	if(tag.compare("") != 0){
	        log_stream << dir << base << "_" << tag << "_log.txt";
	}
	else{
	        log_stream << dir << base << "_log.txt";
	}
        return log_stream.str();
}

string get_edat_name(string base, string dir, string tag){
        ostringstream edat_stream;

	if(tag.compare("") != 0){
	        edat_stream << dir << base << "_" << tag << "_edata.txt";
	}
	else{
	        edat_stream << dir << base << "_edata.txt";
	}
        return edat_stream.str();
}

string get_disp_name(string base, string dir, string tag){
        ostringstream disp_stream;

	if(tag.compare("") != 0){
	        disp_stream << dir << base << "_" << tag;
	}
	else{
	        disp_stream << dir << base;
	}
        return disp_stream.str();
}

vector<string> split(string input, char delim){
    vector<string> result;
    size_t start, iter, len;
    bool reading_token = false;

    for(iter = 0; iter < input.size(); iter++){
        if(delim != input[iter]){
            if(!reading_token){
                reading_token = true;
                start = iter;
            }
        }

        else{
            if(reading_token){
                reading_token = false;
                result.push_back(input.substr(start, iter - start));
            }
        }
    }

    if(reading_token) result.push_back(input.substr(start, iter - start));


    return result;
}

bool get_file_lists(string input, string tag, vector<string> &dat_files, vector<string> &log_files, vector<string> &disp_files, vector<string> &edat_files){

	string dir_name, base_name, name;
	path p;
	vector<string> tokens;
	ostringstream oss, pattern;
	boost::regex expr;
	boost::smatch result;
	int iter;

        if(input.compare("") == 0) return false;

	//Split the response into tokens, delimited by a "/"
	tokens = split(input, '/');
        base_name = *tokens.rbegin();

	if(tokens.size() > 1){
		if(input[0] == '/') oss << "/";
		for(iter = 0; iter < tokens.size() - 1; iter ++){
	        	oss << tokens[iter] << "/";
		}
		dir_name = oss.str();
	}
	else{
		dir_name = "./";
	}
	p = path(dir_name);

	//If a wild card (*) was entered, set the base name to an empty string
	if(base_name.compare("*") == 0) base_name = "";

        //Create the regular expression
	pattern << "^(" << base_name << ".*)" << "\\.dat";
	expr = boost::regex{pattern.str()};

	if(exists(p)){
		for(auto i=directory_iterator(p); i!=directory_iterator(); i++){
			if(! is_directory(i->path())){
				name = i->path().filename().string();
				if(boost::regex_search(name, result, expr)){
					dat_files.push_back(get_dat_name(result[0], dir_name));
					log_files.push_back(get_log_name(result[1], dir_name, tag));
					disp_files.push_back(get_disp_name(result[1], dir_name, tag));
					edat_files.push_back(get_edat_name(result[1], dir_name, tag));

				}
			}
		}
	}

	//Sort the results
        sort(dat_files.begin(), dat_files.end());
        sort(log_files.begin(), log_files.end());
        sort(disp_files.begin(), disp_files.end());
        sort(edat_files.begin(), edat_files.end());


	return true;
}
