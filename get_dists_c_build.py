# file "get_dists_c_build.py"

from cffi import FFI
ffibuilder = FFI()

ffibuilder.cdef("""void pts2pts(float *, const float *, const float *, const int, const int);
                   void cdf(float *, float *, const int, const int);
                """) 

ffibuilder.set_source("_get_dists",
r"""
    #include <math.h>
    #include <float.h>
    
    void pts2pts(float *out, const float *pts_from, const float *pts_to, const int n_from, const int n_to) {
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
                const float dist = (diff_x * diff_x) + (diff_y * diff_y);
                if (dist < p_min) p_min = dist;
            }
            out[i] = sqrtf(p_min);
        }
    }

    static int compare (const void * a, const void * b) {
        float fa = *(const float*) a;
        float fb = *(const float*) b;
        return (fa > fb) - (fa < fb);
    }

    void cdf(float *out, float *ft_dists, const int n_ptsin, const int max_r) {
        const float npts_div = (float) 1.0f / n_ptsin;
        int i, curr;

        out[0] = 0.0f;
        qsort(ft_dists, n_ptsin, sizeof(float), &compare);
        for (i=0, curr=0; i < max_r-1; ++i) { 
            while (ft_dists[curr] < i) {
                out[i] += 1.0f;
                ++curr;
                if (curr == n_ptsin) {
                    int j;
                    for (j=i; j <= max_r; ++j) { 
                        out[j] = n_ptsin * npts_div;
                    }
                    return;
                }
            }
            out[i+1] = out[i];
            out[i] = out[i] * npts_div;
        }
        out[max_r] = out[max_r] * npts_div;
        return;
    }
""")

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)