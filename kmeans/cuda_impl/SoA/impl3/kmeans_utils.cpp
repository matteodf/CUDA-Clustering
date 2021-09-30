
#include "kmeans_utils.h"

void mallocPS(POINT* ps, int n) {
    ps->x = (float*)malloc(sizeof(float) * n);
    ps->y = (float*)malloc(sizeof(float) * n);
    ps->group = (int*)malloc(sizeof(int) * n);
}

void freePS(POINT* ps) {
    free(ps->x);
    free(ps->y);
    free(ps->group);
}

void printDebug(POINT centroids, int num_clusters) {
    printf("centroids: -----------------\n");
    for (int i = 0; i < num_clusters; i++) {
        printf("\t%f %f, %d\n", centroids.x[i], centroids.y[i], centroids.group[i]);
    }
}

void gen_xy(POINT pts, int num_pts, float radius) {
    int i;
    float ang, r;

    for (i = 0; i < num_pts; i++) {
        ang = 2.0 * M_PI * rand() / (RAND_MAX - 1.);
        r = radius * rand() / (RAND_MAX - 1.);
        pts.x[i] = r * cos(ang);
        pts.y[i] = r * sin(ang);
    }
}

int nearest(POINT pts, int ptPos, POINT centroids, int n_cluster) {
    int i, clusterIndex;
    float d, min_d;

    min_d = HUGE_VAL;
    clusterIndex = pts.group[ptPos];
    for (i = 0; i < n_cluster; i++) {
        d = dist2(centroids.x[i], centroids.y[i], pts.x[ptPos], pts.y[ptPos]);
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

void initClusters(POINT pts, int num_pts, POINT centroids, int num_clusters) {
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
    centroids.x[0] = pts.x[selectedIndex];
    centroids.y[0] = pts.y[selectedIndex];
    centroids.group[0] = pts.group[selectedIndex];

    for (j = 0; j < num_pts; ++j)
        shortestDistance[j] = HUGE_VAL;

    /* Select the centroids for the remaining clusters. */
    for (cluster = 1; cluster < num_clusters; cluster++) {
        /* For each point find its closest distance to any of
       the previous cluster centers */
        for (j = 0; j < num_pts; j++) {
            d = dist2(pts.x[j], pts.y[j], centroids.x[cluster - 1], centroids.y[cluster - 1]);

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
        centroids.x[cluster] = pts.x[selectedIndex];
        centroids.y[cluster] = pts.y[selectedIndex];
        centroids.group[cluster] = pts.group[selectedIndex];
    }

    /* Assign each point the index of it's nearest cluster centroid. */
    for (j = 0; j < num_pts; j++)
        pts.group[j] = nearest(pts, j, centroids, num_clusters);

    free(shortestDistance);
    free(cumulativeDistances);
} /* end, initClusters */

void printResult(POINT pts, int num_pts, POINT centroids, int num_clusters) {
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
        if (max_x < pts.x[j])
            max_x = pts.x[j];
        if (min_x > pts.x[j])
            min_x = pts.x[j];
        if (max_y < pts.y[j])
            max_y = pts.y[j];
        if (min_y > pts.y[j])
            min_y = pts.y[j];
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
            if (pts.group[j] != i)
                continue;
            printf("%.3f %.3f c\n", (pts.x[j] - cx) * scale + W / 2,
                   (pts.y[j] - cy) * scale + H / 2);
        }
        printf("\n0 setgray %g %g s\n", (centroids.x[i] - cx) * scale + W / 2,
               (centroids.y[i] - cy) * scale + H / 2);
    }
    printf("\n%%%%EOF");

    free(colors);

    return;
} /* end printResult */