import numpy as np
from numba import guvectorize

# fast matrix multiplication
@guvectorize(['void(float64[:,:], float32[:], float64[:,:])'],
            "(m,n),(n)->(m,n)", nopython=True, cache=True)
def gu_mult(pt, ls, dt):
    for i in range(pt.shape[0]):
        for j in range(pt.shape[1]):
            dt[i,j] = pt[i,j] * ls[j]

# # takes vectors of points and medial axis line segments, returns min distances
# # points = np.array of points (points are lists of length 2)
# # lines  = np.array of pairs of points
# # NOTE: dist_to_lines takes and returns Numpy arrays, not lists.
def dist_pts2lines(points,lines):
    starts = lines[:,0]
    ends = lines[:,1]

    # create a mesh grid to compare all points with all lines 
    point_segs = points[:,np.newaxis,:] - starts[np.newaxis,:,:]

    # calculate dot products [(x1,y1).(x2,y2) = (x1x2+y1y2)]
    # between each point and each line segment (start point)
    line_segs =  ends - starts
    ps_x = np.ascontiguousarray(point_segs[:,:,0], dtype=np.float64)
    ls_x = np.ascontiguousarray(line_segs[:,0], dtype=np.float32)
    ps_y = np.ascontiguousarray(point_segs[:,:,1], dtype=np.float64)
    ls_y = np.ascontiguousarray(line_segs[:,1], dtype=np.float32)

    dots_x = np.empty_like(ps_x)
    dots_y = np.empty_like(dots_x)

    gu_mult(ps_x, ls_x, dots_x)
    gu_mult(ps_y, ls_y, dots_y)
    dots = dots_x + dots_y
    
    # multiply dots by normal vectors to each line to find
    # the projections of each point onto each line
    norms = np.linalg.norm(line_segs, axis=1)
    # lines with norm of 0 represent points and should have projections of 0
    dn = np.divide(dots,norms, out=np.zeros_like(dots), where=norms!=0)
    projections = line_segs * dn[:,:,np.newaxis]

    # find distance from each point to its projection (nearest point on line)
    distances = np.linalg.norm(point_segs-projections, axis=2)
    return distances.min(axis=1) # min distance from each point to all lines




# # takes vectors of points and medial axis line segments, returns min distances
# # points = np.array of points (points are lists of length 2)
# # lines  = np.array of pairs of points
# # NOTE: dist_to_lines takes and returns Numpy arrays, not lists.
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