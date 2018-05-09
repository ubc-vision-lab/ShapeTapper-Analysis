import numpy as np
from numba import guvectorize

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

# Returns distances between 2D sets of points
def points2points(points_from,points_to):
    dists = np.empty(points_from.shape[0], dtype=np.float32)
    gu_dist_2D(points_from, points_to, dists)
    return dists


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