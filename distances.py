import numpy as np
from numba import jit, guvectorize

# Returns distances between 2D sets of points
# If a and b are larger than 2D (p>2), only the first two dims are used
@guvectorize(['void(float32[:,:], float32[:,:], float32[:])'],
              "(m,p),(n,p)->(m)", nopython=True, cache=True)
def gu_dist_2D(a,b,ds):
    for i in range(a.shape[0]) :
        min_d = 1.7976931348623157e+308 # max float
        for j in range(b.shape[0]) :
            x = a[i,0] - b[j,0]
            y = a[i,1] - b[j,1]
            dist = (x*x)+(y*y)
            if (dist < min_d) : min_d = dist
        ds[i] = np.sqrt(min_d)


# Cumulative Density Function for Spatial Analysis
# Takes sorted distances for a set of points (ft) & a vector of distances (rs)
# Counts number of points within each distance (rs[r]), as a function cdf
@guvectorize(['void(float32[:], float32[:], float32[:])'],
              "(m),(n)->(m)", nopython=True, cache=True)
def gu_cdf(rs,ft,cdf):
    idx = 0
    cdf[0] = 0
    npts = ft.shape[0]
    for r in range(rs.shape[0]-1) :
        while ft[idx] < rs[r] :
            cdf[r] += 1
            idx += 1
            if idx == npts :
                cdf[r+1:] = npts
                return
        cdf[r+1] = cdf[r]


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

# Generates n random coordinates within a specified circle
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


# Returns distances between 2D sets of points
def points2points(points_from,points_to):
    pts_from = np.ascontiguousarray(points_from, dtype=np.float32)
    pts_to = np.ascontiguousarray(points_to, dtype=np.float32)
    dists = np.empty(pts_from.shape[0], dtype=np.float32)
    gu_dist_2D(pts_from, pts_to, dists)
    return dists


# cumulative distribution function for "reference object" defined by ro_pts
def cdf(points_in, ro_pts, regions):
    ft_dists = np.sort( points2points(points_in,ro_pts) )
    cdf = np.empty(regions.shape[0], dtype=np.float32) # output
    gu_cdf(regions, ft_dists, cdf)
    return np.true_divide(cdf, points_in.shape[0]) # normalize to [0,1]


# fast matrix-vector multiplication
@guvectorize(['void(float32[:,:], float32[:], float32[:,:])'],
            "(m,n),(n)->(m,n)", nopython=True, cache=True)
def gu_mult_diag(pt, ls, dt):
    for i in range(pt.shape[0]):
        for j in range(pt.shape[1]):
            dt[i,j] = pt[i,j] * ls[j]


# takes vectors of points and medial axis line segments, returns min distances
# points = np.array of points (n x 2)
# lines  = np.array of pairs of points (n x 2 x 2)
def pts2lines(points, lines):
    starts = lines[:,0]
    ends = lines[:,1]
    
    # create a mesh grid to compare all points with all lines 
    point_segs = points[:,np.newaxis,:] - starts[np.newaxis,:,:]
    
    # calculate dot products [(x1,y1).(x2,y2) = (x1x2+y1y2)]
    # between each point and each line segment (start point)
    line_segs =  ends - starts
    ps_x = np.ascontiguousarray(point_segs[:,:,0], dtype=np.float32)
    ls_x = np.ascontiguousarray(line_segs[:,0], dtype=np.float32)
    ps_y = np.ascontiguousarray(point_segs[:,:,1], dtype=np.float32)
    ls_y = np.ascontiguousarray(line_segs[:,1], dtype=np.float32)
    dots_x = np.empty_like(ps_x)
    dots_y = np.empty_like(dots_x)

    gu_mult_diag(ps_x, ls_x, dots_x)
    gu_mult_diag(ps_y, ls_y, dots_y)
    dots = dots_x + dots_y
    
    # avoid taking the square root when finding norms
    norms_sq = np.sum(np.square(line_segs), axis=1)
    
    # lines with norm of 0 represent points and should have projections of 0
    # constrain dot products b/t 0 and 1 to handle projections outside of line
    dn = np.divide(dots, norms_sq, out=np.zeros_like(dots), where=norms_sq!=0)
    np.clip(dn, 0., 1., out=dn)
    
    projections = line_segs * dn[:,:,np.newaxis]

    # find distance from each point to its projection (nearest point on line)
    distances = np.linalg.norm(point_segs-projections, axis=2)
    return distances # min distance from each point to all lines