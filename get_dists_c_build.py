# file "get_dists_c_build.py"

from cffi import FFI
ffibuilder = FFI()

ffibuilder.cdef("""void pts2pts(float *, const float *, const float *, int, int);
                   void cdf(int *, float *, int, const float *, int);
                   void gen_uniform_bounds(float *, const float *, int, float, float, float, float, int, int);
                   void gen_uniform_circle(float *, float, float, float, int, int);
                """) 

ffibuilder.set_source("_get_dists",
r"""    
    #include <math.h>
    #include <float.h>
    #include "spatlib.c"
""")

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)