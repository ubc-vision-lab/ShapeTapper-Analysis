import numpy as np
from numba import guvectorize


@guvectorize(['void(float64[:], float64[:], float64[:,:])'],
              "(m),(n)->(m,n)", nopython=True, cache=True)
def gu_mult_dist(a,b,ds):
    for i in range(a.shape[0]):
        for j in range(b.shape[0]):
            ds[i,j] = np.square(a[i] - b[j])
    
    
@guvectorize(['void(float64[:,:], float64[:,:], float64[:,:])'],
              "(m,n),(m,n)->(m,n)", nopython=True, cache=True)   
def gu_sum(a, b, s):
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            s[i,j] = a[i,j] + b[i,j]


@guvectorize(['void(float64[:,:], float64[:,:])'],
              "(m,n)->(m,n)", nopython=True, cache=True)
def gu_sqrt(x, sqt):
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            sqt[i,j] = np.sqrt(x[i,j])


# fast matrix-vector multiplication
@guvectorize(['void(float64[:,:], float64[:], float64[:,:])'],
            "(m,n),(n)->(m,n)", nopython=True, cache=True)
def gu_mult_diag(pt, ls, dt):
    for i in range(pt.shape[0]):
        for j in range(pt.shape[1]):
            dt[i,j] = pt[i,j] * ls[j]



# takes vector of points and centroid, returns distances from centroid
# points = np.array of points (points are lists of length 2)
# centroid = list of centroid coordinates
def points2points(points_from,points_to):
    from_y = np.ascontiguousarray(points_from[:,0], dtype=np.float64)
    from_x = np.ascontiguousarray(points_from[:,1], dtype=np.float64)
    to_y = np.ascontiguousarray(points_to[:,0], dtype=np.float64)
    to_x = np.ascontiguousarray(points_to[:,1], dtype=np.float64)
    
    ds_y = np.empty((from_y.size,to_y.size), dtype=np.float64)
    ds_x = np.empty((from_x.size,to_x.size), dtype=np.float64)
    ds_sum = np.empty((from_x.size,to_x.size), dtype=np.float64)
    dists = np.empty((from_x.size,to_x.size), dtype=np.float64)
    
    gu_mult_dist(from_y,to_y,ds_y)
    gu_mult_dist(from_x,to_x,ds_x)
    gu_sum(ds_x, ds_y, ds_sum)
    gu_sqrt(ds_sum, dists)

    return dists
#    point_diffs = points_from[:,np.newaxis,:] - points_to[np.newaxis,:,:]
#    return np.linalg.norm(point_diffs, axis=2)


# takes vectors of points and medial axis line segments, returns min distances
# points = np.array of points (points are lists of length 2)
# lines  = np.array of pairs of points
def pts2lines(points,lines):
    starts = lines[:,0]
    ends = lines[:,1]
    
    # create a mesh grid to compare all points with all lines 
    point_segs = points[:,np.newaxis,:] - starts[np.newaxis,:,:]
    
    # calculate dot products [(x1,y1).(x2,y2) = (x1x2+y1y2)]
    # between each point and each line segment (start point)
    line_segs =  ends - starts
    
    ps_x = np.ascontiguousarray(point_segs[:,:,0], dtype=np.float64)
    ls_x = np.ascontiguousarray(line_segs[:,0], dtype=np.float64)
    ps_y = np.ascontiguousarray(point_segs[:,:,1], dtype=np.float64)
    ls_y = np.ascontiguousarray(line_segs[:,1], dtype=np.float64)
    
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
    np.clip(dn,0.,1., out=dn)
    
    projections = line_segs * dn[:,:,np.newaxis]

    # find distance from each point to its projection (nearest point on line)
    distances = np.linalg.norm(point_segs-projections, axis=2)
    return distances # min distance from each point to all lines


# takes vector of points and centroid, returns distances from centroid
# points = np.array of points (points are lists of length 2)
# centroid = list of centroid coordinates
def points2point(points,target):
    return np.linalg.norm(points - target[0:2],axis=1)


## takes vectors of points and medial axis line segments, returns min distances
## points = np.array of points (points are lists of length 2)
## lines  = np.array of pairs of points
## NOTE: dist_to_lines takes and returns Numpy arrays, not lists.
# def dist_pts2lines(points,lines):
#     starts = lines[:,0]
#     ends = lines[:,1]

#     # create a mesh grid to compare all points with all lines 
#     point_segs = points[:,np.newaxis,:] - starts[np.newaxis,:,:]

#     # calculate dot products [(x1,y1).(x2,y2) = (x1x2+y1y2)]
#     # between each point and each line segment (start point)
#     line_segs =  ends - starts
#     line_segs_diag_x = np.diag(line_segs[:,0])
#     line_segs_diag_y = np.diag(line_segs[:,1])
#     dots = ( np.matmul(point_segs[:,:,0],line_segs_diag_x) + 
#                  np.matmul(point_segs[:,:,1],line_segs_diag_y) ) 
    
#     # multiply dots by normal vectors to each line to find
#     # the projections of each point onto each line
#     norms = np.linalg.norm(line_segs, axis=1)
#     dn = dots/norms
#     projections = line_segs * dn[:,:,np.newaxis]
    
#     # find distance from each point to its projection (nearest point on line)
#     distances = np.linalg.norm(point_segs-projections, axis=2)
#     return distances.min(axis=1) # min distance from each point to all lines


# # same as dis_pts2lines, but takes a 3D matrix of points which represent 
# # simulated sets of point observations, calculations must be done in 4-D
# def dist_pts2lines_4d(points,lines):
#     starts = lines[:,0]
#     ends = lines[:,1]

#     # starts must be tiled to 3-D to create a 4-D mesh grid
#     dummy = np.zeros(points[0].shape)
#     starts_t = dummy[:,np.newaxis,:] + starts[np.newaxis,:,:]
#     point_segs = points[:,:,np.newaxis,:] - starts_t[np.newaxis,:,:,:]

#     # calculate dot products [(x1,y1).(x2,y2) = (x1x2+y1y2)]
#     # between each point and each line segment (start point)
#     line_segs =  ends - starts
#     ls_diag_x = np.diag(line_segs[:,0])
#     ls_diag_y = np.diag(line_segs[:,1])
#     dots = ( np.matmul(point_segs[:,:,:,0],ls_diag_x[np.newaxis,:,:]) + 
#                  np.matmul(point_segs[:,:,:,1],ls_diag_y[np.newaxis,:,:]) )
    
#     # multiply dots by norms to each line to find projections
#     norms = np.linalg.norm(line_segs, axis=1)
#     dn = dots / norms
#     lines_t = dummy[:,np.newaxis,:] + line_segs[np.newaxis,:,:]
#     projections = lines_t[np.newaxis,:,:,:] * dn[:,:,:,np.newaxis]
    
#     # find distance from point to each projection
#     distances = np.linalg.norm(point_segs-projections, axis=3)
#     return distances.min(axis=2) # min distance from each point to all lines