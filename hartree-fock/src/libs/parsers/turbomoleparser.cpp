#include "turbomoleparser.h"

#include <fstream>
#include <iostream>
#include <boost/regex.hpp>

using namespace std;
using namespace boost;

TurboMoleParser::TurboMoleParser()
{
}

void TurboMoleParser::read(string fileName)
{
    ifstream dataFile(fileName);
    string line;
    try {
        while (getline(dataFile, line))
        {
            bool skip = false;
            skip |= regex_match(line, regex("#.*"));
            skip |= regex_match(line, regex("$basis.*"));
            skip |= regex_match(line, regex("$end.*"));
            skip |= regex_match(line, regex("\\*.*"));
            if(skip) {
                continue;
            }
            cout << line;
            if(regex_match(line, regex("\\s*n\\s*\\d-\\d+[G]\\s*"))) { // n 4-31G
                cout << "               <--  basis set type";
            }
            if(regex_match(line, regex("\\s*[0-9]\\s*[spdf]\\s*"))) { // 4 s
                cout << "               <--  orbital type";
            }
            if(regex_match(line, regex("\\s*\\d+\\.?\\d+\\s*\\d+\\.?\\d+\\s*"))) { // 2.1 4.9
                cout << "               <--  expontent weight";
            }
            cout << endl;
        }
    } catch (regex_error& error) {
        cout << "Regex error " << error.code() << " " << error.what() << endl;
        cout << (error.code() == regex_constants::error_escape) << endl;
    }
}
