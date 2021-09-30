#ifndef INPUTREADER_H
#define INPUTREADER_H

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Points.h"
#include "utils.h"

class InputReader {
   private:
    std::ifstream fin;
    Points points;
    int n;
    int findNumRows(std::string filename);

   public:
    InputReader(std::string filename);
    void parse();
    Points getPoints();
};

#endif