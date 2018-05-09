import numpy as np
from numba import jit

# Returns -1 if is right of, 1 if is left of
@jit(["float32(float32,float32,float32,float32,float32,float32)"], nopython=True, cache=True)
def is_left_of(v1_x, v1_y, v2_x, v2_y, p_x, p_y) :
    return ( (v2_x-v1_x) * (p_y-v1_y) - (p_x-v1_x) * (v2_y-v1_y) )

# called by gen_uniform_points_bounds()
# performs test using winding number algorithm (no trig functions)
@jit(nopython=True, cache=True)
def point_in_poly(p_x, p_y, edge_pts) :
    wd_num = 0
    for i in range(1,edge_pts.shape[0]) :
        v1_x = edge_pts[i-1,0]
        v1_y = edge_pts[i-1,1]
        v2_x = edge_pts[i,0]
        v2_y = edge_pts[i,1]
        if ( v1_y <= p_y ) :
            if ( v2_y > p_y ) :    # if v1-v2 is upward crossing of ray from p
                if (is_left_of(v1_x, v1_y, v2_x, v2_y, p_x, p_y) > 0) :
                    wd_num += 1
        else :
            if ( v2_y <= p_y ) :   # if v1-v2 is downward crossing of ray from p
                if (is_left_of(v1_x, v1_y, v2_x, v2_y, p_x, p_y) < 0) :
                    wd_num -= 1
    return wd_num

# Generates n random coordinates within a specified polygon given by edge_points
@jit(nopython=True, cache=True)
def gen_uniform_bounds(edge_points, min_x, max_x, min_y, max_y, n_pts, rands) :
    for i in range(n_pts) :
        while (1) :
            rx = np.random.uniform(min_x, max_x)
            ry = np.random.uniform(min_y, max_y)
            if point_in_poly(rx, ry, edge_points) :
                rands[i,0] = rx
                rands[i,1] = ry
                break


# Generates n random coordinates within a specified circle
# Circle: center (x,y) = (dims[0], dims[1]) , radius = dims[2]
@jit(nopython=True, cache=True)
def gen_uniform_circle(cent_x, cent_y, rad, n_pts, rands) :
    min_x  = cent_x - rad
    min_y  = cent_y - rad
    max_x  = cent_x + rad
    max_y  = cent_y + rad
    rad_sq = rad * rad
    
    for i in range(n_pts) :
        while (1) :
            rx = np.random.uniform(min_x, max_x)
            ry = np.random.uniform(min_y, max_y)
            dx = rx - cent_x
            dy = ry - cent_y
            if not ( (dx*dx) + (dy*dy) > rad_sq ) :
                rands[i,0] = rx
                rands[i,1] = ry
                break