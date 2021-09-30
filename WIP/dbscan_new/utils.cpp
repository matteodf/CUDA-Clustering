#include "utils.h"

void printResult(Points pts, int num_pts, int num_clusters) {
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
            if (pts.cluster[j] != i)
                continue;
            printf("%.3f %.3f c\n", (pts.x[j] - cx) * scale + W / 2,
                   (pts.y[j] - cy) * scale + H / 2);
        }
    }

    printf("%g %g %g setrgbcolor\n", 0.85, 0.85,
           0.85);
    for (i = 0; i < num_pts; i++) {
        if (pts.cluster[i] == NOISE) {
            printf("%.3f %.3f c\n", (pts.x[i] - cx) * scale + W / 2,
                   (pts.y[i] - cy) * scale + H / 2);
        }
    }
    printf("\n%%%%EOF");

    free(colors);

    return;
} /* end printResult */