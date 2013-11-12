#ifndef TURBOMOLEPARSER_H
#define TURBOMOLEPARSER_H

#include <string>

using namespace std;

class TurboMoleParser
{
public:
    TurboMoleParser();

    void read(string fileName);
};

#endif // TURBOMOLEPARSER_H
