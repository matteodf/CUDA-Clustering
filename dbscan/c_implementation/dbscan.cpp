/**
 * Sequential implementation of DBScan clustering method
 * 
 * @file dbscan.cpp
 * @author Matteo De Filippis / Michele Andreata
 * @version 0.1
 * 
 */

#include <sys/time.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <sys/stat.h>
//#include <sys/types.h>
//#include <algorithm>
//#include <map>
//#include <string>
//include <bits/stdc++.h>

using namespace std;

const int NOISE = -2;
const int NOT_CLASSIFIED = -1;

/**
    Function used to determine current timestamps in seconds
    @return A double value that indicates the current timestamp in seconds
*/
double seconds() {
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

class Point {
   public:
    double x, y;
    int ptsCnt, cluster;
    double getDis(const Point& ot) {
        return sqrt((x - ot.x) * (x - ot.x) + (y - ot.y) * (y - ot.y));
    }
};

class DBSCAN {
   public:
    int minPts;
    double eps;
    vector<Point> points;
    int size;
    vector<vector<int> > adjPoints;
    vector<bool> visited;
    vector<vector<int> > cluster;
    int clusterIdx;

    DBSCAN(double eps, int minPts, vector<Point> points) {
        this->eps = eps;
        this->minPts = minPts;
        this->points = points;
        this->size = (int)points.size();
        adjPoints.resize(size);
        this->clusterIdx = -1;
    }
    void run() {
        checkNearPoints();

        for (int i = 0; i < size; i++) {
            if (points[i].cluster != NOT_CLASSIFIED) continue;

            if (isCoreObject(i)) {
                dfs(i, ++clusterIdx);
            } else {
                points[i].cluster = NOISE;
            }
        }

        cluster.resize(clusterIdx + 1);
        for (int i = 0; i < size; i++) {
            if (points[i].cluster != NOISE) {
                cluster[points[i].cluster].push_back(i);
            }
        }
    }

    void dfs(int now, int c) {
        points[now].cluster = c;
        if (!isCoreObject(now)) return;

        for (auto& next : adjPoints[now]) {
            if (points[next].cluster != NOT_CLASSIFIED && points[next].cluster != NOISE) continue;
            dfs(next, c);
        }
    }

    void checkNearPoints() {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i == j) continue;
                if (points[i].getDis(points[j]) <= eps) {
                    points[i].ptsCnt++;
                    adjPoints[i].push_back(j);
                }
            }
        }
    }
    // is idx'th point core object?
    bool isCoreObject(int idx) {
        return points[idx].ptsCnt >= minPts;
    }

    vector<vector<int> > getCluster() {
        return cluster;
    }
};

class InputReader {
   private:
    ifstream fin;
    vector<Point> points;

   public:
    InputReader(string filename) {
        filename = "../../input_data/" + filename;
        fin.open(filename);
        if (!fin) {
            cout << filename << " file could not be opened\n";
            exit(0);
        }
        parse();
    }
    void parse() {
        int idx;
        double x, y;
        while (!fin.eof()) {
            fin >> idx >> x >> y;
            points.push_back({x, y, 0, NOT_CLASSIFIED});
        }
        points.pop_back();
    }
    vector<Point> getPoints() {
        return points;
    }
};


void printResults(DBSCAN dbScan, string filename) {
            #define W 400
            #define H 400

        
            size_t lastdot = filename.find_last_of(".");
            if (lastdot != std::string::npos)
                filename = filename.substr(0, lastdot); 
            
            FILE * pFile;

            
            // Creating a directory
            if (mkdir("../../output_data/dbScan_c_impl", 0777) != -1)
                cout << "New directory created for DBscan - C Implementation" << endl;

            string outputFilename = "../../output_data/dbScan_c_impl/" + filename + ".eps";

            pFile = fopen(outputFilename.c_str(), "w");


            // remove ".txt" from filename
            if(filename.size()<4){
                cout << filename << "input file name's format is wrong\n";
                exit(0);
            }
            for(int i=0;i<4;i++) 
                filename.pop_back();
            

            int i, j;
            double min_x, max_x, min_y, max_y, scale, cx, cy;
            double* colors = (double*)malloc(sizeof(double) * (dbScan.clusterIdx + 1) * 3);

            for (i = 0; i < dbScan.clusterIdx + 1; i++) {
                colors[3 * i + 0] = (3 * (i + 1) % 11) / 11.;
                colors[3 * i + 1] = (7 * i % 11) / 11.;
                colors[3 * i + 2] = (9 * i % 11) / 11.;
            }

            max_x = max_y = -HUGE_VAL;
            min_x = min_y = HUGE_VAL;
            for (j = 0; j < dbScan.size; j++) {
                if (max_x < dbScan.points[j].x)
                    max_x = dbScan.points[j].x;
                if (min_x > dbScan.points[j].x)
                    min_x = dbScan.points[j].x;
                if (max_y < dbScan.points[j].y)
                    max_y = dbScan.points[j].y;
                if (min_y > dbScan.points[j].y)
                    min_y = dbScan.points[j].y;
            }

            scale = W / (max_x - min_x);
            if (scale > H / (max_y - min_y))
                scale = H / (max_y - min_y);
            cx = (max_x + min_x) / 2;
            cy = (max_y + min_y) / 2;

            fprintf(pFile, "%%!PS-Adobe-3.0\n%%%%BoundingBox: -5 -5 %d %d\n", W + 10, H + 10);
            fprintf(pFile, 
                "/l {rlineto} def /m {rmoveto} def\n"
                "/c { .25 sub exch .25 sub exch .5 0 360 arc fill } def\n"
                "/s { moveto -2 0 m 2 2 l 2 -2 l -2 -2 l closepath "
                "	gsave 1 setgray fill grestore gsave 3 setlinewidth"
                " 1 setgray stroke grestore 0 setgray stroke }def\n");

            for (i = 0; i < dbScan.clusterIdx + 1; i++) {
                fprintf(pFile, "%g %g %g setrgbcolor\n", colors[3 * i], colors[3 * i + 1],
                    colors[3 * i + 2]);

                for (j = 0; j < dbScan.size; j++) {
                    if (dbScan.points[j].cluster != i)
                        continue;
                    fprintf(pFile, "%.3f %.3f c\n", (dbScan.points[j].x - cx) * scale + W / 2,
                        (dbScan.points[j].y - cy) * scale + H / 2);
                }
            }

            fprintf(pFile, "%g %g %g setrgbcolor\n", 0.85, 0.85,
                0.85);
            for (i = 0; i < dbScan.size; i++) {
                if (dbScan.points[i].cluster == NOISE) {
                    fprintf(pFile, "%.3f %.3f c\n", (dbScan.points[i].x - cx) * scale + W / 2,
                        (dbScan.points[i].y - cy) * scale + H / 2);
                }
            }
            fprintf(pFile, "\n%%%%EOF");

            free(colors);

            return;
        } /* end printResult */


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
