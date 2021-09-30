#include "kmeans_utils.h"

#include <fstream>

void printDebug(POINT* centroids, int num_clusters) {
    printf("centroids: -----------------\n");
    for (int i = 0; i < num_clusters; i++) {
        printf("\t%f %f, %d\n", centroids[i].x, centroids[i].y, centroids[i].group);
    }
}

void genRandomPoints(POINT* pts, int num_pts, float radius) {
    int i;
    float ang, r;

    for (i = 0; i < num_pts; i++) {
        ang = 2.0 * M_PI * rand() / (RAND_MAX - 1.);
        r = radius * rand() / (RAND_MAX - 1.);
        pts[i].x = r * cos(ang);
        pts[i].y = r * sin(ang);
    }
}

int findNumRows(char* filename) {
    int count = 0;
    std::ifstream in(filename);
    std::string unused;
    while (std::getline(in, unused))
        ++count;
    in.close();
    return count;
}

void readPointsFromFile(POINT* pts, char* filename) {
    std::ifstream fin;
    fin.open(filename);
    if (!fin) {
        printf("%s file could not be opened\n", filename);
        exit(0);
    }
    int idx;
    float x, y;
    while (!fin.eof()) {
        fin >> idx >> x >> y;
        pts[idx].x = x;
        pts[idx].y = y;
    }
    fin.close();
}

int nearest(POINT* pt, POINT* cent, int n_cluster) {
    int i, clusterIndex;
    float d, min_d;

    min_d = HUGE_VAL;
    clusterIndex = pt->group;
    for (i = 0; i < n_cluster; i++) {
        d = dist2(&cent[i], pt);
        if (d < min_d) {
            min_d = d;
            clusterIndex = i;
        }
    }
    return clusterIndex;
}

int bisectionSearch(float* x, int n, float v) {
    int il, ir, i;

    if (n < 1) {
        return 0;
    }
    /* If v is less than x(0) or greater than x(n-1)  */
    if (v < x[0]) {
        return 0;
    } else if (v > x[n - 1]) {
        return n - 1;
    }

    /*bisection search */
    il = 0;
    ir = n - 1;

    i = (il + ir) / 2;
    while (i != il) {
        if (x[i] <= v) {
            il = i;
        } else {
            ir = i;
        }
        i = (il + ir) / 2;
    }

    if (x[i] <= v)
        i = ir;
    return i;
}

void initClustersRandom(POINT* pts, int num_pts, int num_clusters) {
    for (int i = 0; i < num_pts; i++) {
        pts[i].group = i % num_clusters;
    }
}

void initClusters(POINT* pts, int num_pts, POINT* centroids, int num_clusters) {
    int j;
    int selectedIndex;
    int cluster;
    float sum;
    float d;
    float random;
    float* cumulativeDistances;
    float* shortestDistance;

    cumulativeDistances = (float*)malloc(sizeof(float) * num_pts);
    shortestDistance = (float*)malloc(sizeof(float) * num_pts);

    /* Pick the first cluster centroids at random. */
    selectedIndex = rand() % num_pts;
    centroids[0] = pts[selectedIndex];

    for (j = 0; j < num_pts; ++j)
        shortestDistance[j] = HUGE_VAL;

    /* Select the centroids for the remaining clusters. */
    for (cluster = 1; cluster < num_clusters; cluster++) {
        /* For each point find its closest distance to any of
       the previous cluster centers */
        for (j = 0; j < num_pts; j++) {
            d = dist2(&pts[j], &centroids[cluster - 1]);

            if (d < shortestDistance[j])
                shortestDistance[j] = d;
        }

        /* Create an array of the cumulative distances. */
        sum = 0.0;
        for (j = 0; j < num_pts; j++) {
            sum += shortestDistance[j];
            cumulativeDistances[j] = sum;
        }

        /* Select a point at random. Those with greater distances
        have a greater probability of being selected. */
        random = (float)rand() / (float)RAND_MAX * sum;
        selectedIndex = bisectionSearch(cumulativeDistances, num_pts, random);

        /* assign the selected point as the center */
        centroids[cluster] = pts[selectedIndex];
    }

    /* Assign each point the index of it's nearest cluster centroid. */
    for (j = 0; j < num_pts; j++)
        pts[j].group = nearest(&pts[j], centroids, num_clusters);

    free(shortestDistance);
    free(cumulativeDistances);
} /* end, initClusters */

void writePointsToFile(char* filename, POINT* pts, int num_pts, POINT* centroids, int num_clusters) {
    FILE* fp;
    fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("file can't be opened\n");
        exit(1);
    }
    for (int i = 0; i < num_pts; i++) {
        fprintf(fp, "%f %f %d\n", pts[i].x, pts[i].y, pts[i].group);
    }
    for (int i = 0; i < num_clusters; i++) {
        fprintf(fp, "%f %f %d c\n", centroids[i].x, centroids[i].y, centroids[i].group);
    }
    fclose(fp);
}

void printResult(POINT* pts, int num_pts, POINT* centroids, int num_clusters) {
#define W 400
#define H 400

    int i, j;
    double min_x, max_x, min_y, max_y, scale, cx, cy;
    double* colors = (double*)malloc(sizeof(double) * num_clusters * 3);

    for (i = 0; i < num_clusters; i++) {
        colors[3 * i + 0] = (3 * (i + 1) % 11) / 11.;
        colors[3 * i + 1] = (7 * i % 11) / 11.;
        colors[3 * i + 2] = (9 * i % 11) / 11.;
    }

    max_x = max_y = -HUGE_VAL;
    min_x = min_y = HUGE_VAL;
    for (j = 0; j < num_pts; j++) {
        if (max_x < pts[j].x)
            max_x = pts[j].x;
        if (min_x > pts[j].x)
            min_x = pts[j].x;
        if (max_y < pts[j].y)
            max_y = pts[j].y;
        if (min_y > pts[j].y)
            min_y = pts[j].y;
    }

    scale = W / (max_x - min_x);
    if (scale > H / (max_y - min_y))
        scale = H / (max_y - min_y);
    cx = (max_x + min_x) / 2;
    cy = (max_y + min_y) / 2;

    printf("%%!PS-Adobe-3.0\n%%%%BoundingBox: -5 -5 %d %d\n", W + 10, H + 10);
    printf(
        "/l {rlineto} def /m {rmoveto} def\n"
        "/c { .25 sub exch .25 sub exch .5 0 360 arc fill } def\n"
        "/s { moveto -2 0 m 2 2 l 2 -2 l -2 -2 l closepath "
        "	gsave 1 setgray fill grestore gsave 3 setlinewidth"
        " 1 setgray stroke grestore 0 setgray stroke }def\n");

    for (i = 0; i < num_clusters; i++) {
        printf("%g %g %g setrgbcolor\n", colors[3 * i], colors[3 * i + 1],
               colors[3 * i + 2]);

        for (j = 0; j < num_pts; j++) {
            if (pts[j].group != i)
                continue;
            printf("%.3f %.3f c\n", (pts[j].x - cx) * scale + W / 2,
                   (pts[j].y - cy) * scale + H / 2);
        }
        printf("\n0 setgray %g %g s\n", (centroids[i].x - cx) * scale + W / 2,
               (centroids[i].y - cy) * scale + H / 2);
    }
    printf("\n%%%%EOF");

    free(colors);

    return;
} /* end printResult */