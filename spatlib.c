#include <stdlib.h>
#include <math.h>
#include <float.h>

// calculates minimum distance of each points-from to the set of points-to
// NOTE: this function crashes for n_from > 250, so there is a fallback NumPy function in distances.py
void pts2pts(float *out, const float *pts_from, const float *pts_to, int n_from, int n_to) {
    int i,j;
    for (i=0; i<n_from; ++i) 
    {
        float p_min = FLT_MAX;
        const float from_x = pts_from[i*2];
        const float from_y = pts_from[i*2+1]; 
        for (j=0; j<n_to; ++j) 
        {
            const float diff_x = from_x - pts_to[j*2];
            const float diff_y = from_y - pts_to[j*2+1];
            const float dist = (diff_x * diff_x) + (diff_y * diff_y); // (x2-x1)^2 + (y2-y1)^2
            if (dist < p_min) p_min = dist; // avoid using sqrt until after comparisons
        }
        out[i] = sqrtf(p_min);
    }
}


// helper function for qsort in cdf()
// returns -1 for a<b, 0 for equality and 1 for a>b 
static int compare (const void * a, const void * b) {
    const float fa = *(const float*) a;
    const float fb = *(const float*) b;
    return (fa > fb) - (fa < fb);
}

// generates a CDF for feature-transformed set of points (dist from each point to reference object) 
// at a increasing radius (r) from the reference object (ratio: #points in radius to total #points)
void cdf(int *out, float *ft_dists, int n_ptsin, int max_r) {
    int r, curr;
    out[0] = 0;

    qsort(ft_dists, n_ptsin, sizeof(float), &compare);
    
    for (r=0, curr=0; r < max_r-1; r++) { 
        while (ft_dists[curr] < r) {
            out[r]++;
            curr++;
            if (curr == n_ptsin) {
                int i;
                for (i=r; i < max_r; ++i) { 
                    out[i] = n_ptsin;
                }
                return;
            }
        }
        out[r+1] = out[r];
    }
    return;
}


// generates n float in range (min,max)
static float randf(float min, float max) {
    if (min == 0.0f) {
        return (float) ((double)rand() / (double)(RAND_MAX) * (double) (max)) ;
    } else {
        const double min_d = (double) min ;
        const double range = (double) max - min_d;
        return (float) (min_d + ((double)rand() / (double)(RAND_MAX)) * range) ;
    }
}

// helper for in_polygon_test()
static int is_left_of(const float v1_x, const float v1_y, const float v2_x, const float v2_y, float p_x, float p_y) {
    return (int) ( (v2_x - v1_x) * (p_y - v1_y) - (p_x -  v1_x) * (v2_y - v1_y) );
}

// called by gen_uniform_points_bounds()
// performs test using winding number algorithm (no trig functions)
static int point_in_poly(float p_x, float p_y, const float * edge_pts, int n_edges) {
    int wd_num = 0;
    int i;
    float v1_x, v1_y, v2_x, v2_y;
    for (i=1; i<n_edges; ++i) {
        v1_x = edge_pts[(i-1)*2];
        v1_y = edge_pts[(i-1)*2+1];
        v2_x = edge_pts[i*2];
        v2_y = edge_pts[i*2+1];
        if ( v1_y <= p_y ) {
            if ( v2_y > p_y ) {     // if v1-v2 is upward crossing of ray from p
                if (is_left_of(v1_x, v1_y, v2_x, v2_y, p_x, p_y) > 0) {
                    wd_num++;
                }
            }
        } else {
            if ( v2_y <= p_y ) { // if v1-v2 is downward crossing of ray from p
                if (is_left_of(v1_x, v1_y, v2_x, v2_y, p_x, p_y) < 0) {
                    wd_num--;
                }
            }
        }
    }
    return wd_num;
}

// generates n random coordinates within an arbitrary polygon
void gen_uniform_bounds(float * out, const float * edge_points, int n_edges, float min_x, float max_x, float min_y, float max_y, int n_pts, int seed) {
    int num_unif = 0;
    float rx, ry;

    srand(seed); // init random number generator

    while (1) {
        rx = randf(min_x, max_x);
        ry = randf(min_y, max_y);
        if ( point_in_poly(rx, ry, edge_points, n_edges) ) {
            out[num_unif*2]   = rx;
            out[num_unif*2+1] = ry;
            num_unif++;
        }
        if (num_unif == n_pts) {
            return;
        }
    }
}

// generates n random coordinates within a specified circle
void gen_uniform_circle(float * out, float cent_x, float cent_y, float rad, int n_pts, int seed) {
    int num_unif = 0;
    float rx, ry, dx, dy;

    const float rad_sq = rad * rad;
    const float min_x = floorf(cent_x - rad);
    const float min_y = floorf(cent_y - rad);
    const float max_x = ceilf(cent_x + rad);
    const float max_y = ceilf(cent_y + rad);

    srand(seed); // init random number generator

    while (1) {
        rx = randf(min_x, max_x);
        ry = randf(min_y, max_y);
        dx = rx - cent_x;
        dy = ry - cent_y;
        if ( (dx*dx) + (dy*dy) <= rad_sq ) {
            out[num_unif*2]   = rx;
            out[num_unif*2+1] = ry;
            num_unif++;
        }
        if (num_unif == n_pts) {
            return;
        }
    }
}