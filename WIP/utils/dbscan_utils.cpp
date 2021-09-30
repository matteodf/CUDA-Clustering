#include "dbscan_utils.h"

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

/**
    Function used to print results
    @param dbScan a DBSCAN structure containing number of total points and a structure of points defined by x and y coordinates and their clusterIdx
    @param filename filename of the input file, used to create the related output file
*/
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
