#include "InputReader.h"
#include "Points.h"
#include "dbscan.h"
#include "utils.h"

int main(int argc, const char* argv[]) {
    if (argc != 4) {
        std::cout << "Please follow this format. clustering.exe [intput] [eps] [minPts]\n";
        return 0;
    }

    std::string inputFileName(argv[1]);
    std::string eps(argv[2]);
    std::string minPts(argv[3]);

    InputReader inputReader(inputFileName);

    DBSCAN dbScan(stod(eps), stoi(minPts), inputReader.getPoints());
    double start = seconds();
    dbScan.run();
    double stop = seconds();
    if (!PRINT_RESULTS) printf("DBSCAN total time: %f\n", stop - start);
    if (!PRINT_RESULTS) printf("num of clusters: %d\n", dbScan.n);

    if (PRINT_RESULTS) printResult(dbScan.points, dbScan.size, dbScan.n);
    return 0;
}