#include "InputReader.h"

InputReader::InputReader(std::string filename) {
    fin.open(filename);
    if (!fin) {
        std::cout << filename << " file could not be opened\n";
        exit(0);
    }
    n = findNumRows(filename);
    points.allocate(n);
    parse();
}

void InputReader::parse() {
    int idx;
    double x, y;
    while (!fin.eof()) {
        fin >> idx >> x >> y;
        if (idx < n) {
            points.x[idx] = x;
            points.y[idx] = y;
            points.ptsCnt[idx] = 0;
            points.cluster[idx] = NOT_CLASSIFIED;
        }
    }
}
Points InputReader::getPoints() {
    return points;
}

int InputReader::findNumRows(std::string filename) {
    int count = 0;
    std::ifstream in(filename);
    std::string unused;
    while (std::getline(in, unused))
        ++count;
    in.close();
    return count;
}