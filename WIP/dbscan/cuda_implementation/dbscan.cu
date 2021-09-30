/**
 * Parallel implementation of DBScan clustering method
 * 
 * @file dbscan.cu
 * @author Matteo De Filippis / Michele Andreata
 * @version 0.1
 * 
 */

#include "../../utils/dbscan_utils.cu"

int main(int argc, const char* argv[]) {
    if (argc != 4) {
        cout << "Please follow this format. clustering.exe [input] [eps] [minPts]";
        return 0;
    }

    string inputFileName(argv[1]);
    string eps(argv[2]);
    string minPts(argv[3]);

    InputReader inputReader(inputFileName);

    DBSCAN dbScan(stod(eps), stoi(minPts), inputReader.getPoints());
    double start = seconds();
    dbScan.run();
    double stop = seconds();
    printf("DBSCAN time: %f\n", stop - start);
    printf("Number of clusters created: %d\n", dbScan.clusterIdx+1);

    printResults(dbScan, inputFileName);

    return 0;
}
