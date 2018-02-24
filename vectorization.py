import numpy as np

# takes vectors of points and medial axis line segments, returns min distances
# points = np.array of points (points are lists of length 2)
# lines  = np.array of pairs of points
# NOTE: dist_to_lines takes and returns Numpy arrays, not lists.
def dist_pts2lines(points,lines):
    starts = lines[:,0]
    ends = lines[:,1]
    line_segs =  ends - starts

    # points must be repeated in order to compare to each line segment
    pts_exp = np.expand_dims(points,2) #add a new dimension for tiling
    pts_tiled = np.tile(pts_exp, starts.shape[0]).transpose((0,2,1))
    starts_exp = np.expand_dims(starts,2)
    starts_tiled = np.tile(starts_exp, points.shape[0]).transpose((2,0,1))

    point_segs = pts_tiled - starts_tiled

    # calculate dot products [(x1,y1).(x2,y2) = (x1x2+y1y2)]
    line_segs_diag_x = np.diag(line_segs[:,0])
    line_segs_diag_y = np.diag(line_segs[:,1])
    dots = ( np.matmul(point_segs[:,:,0],line_segs_diag_x) + 
                 np.matmul(point_segs[:,:,1],line_segs_diag_y) ) 
    
    norms = np.linalg.norm(line_segs, axis=1)
    dn = np.tile(np.expand_dims(dots/norms,2),2) # tile to match point_segs

    lsegs_exp = np.expand_dims(line_segs,2) # tile to match point_segs
    lsegs_tiled = np.tile(lsegs_exp, points.shape[0]).transpose((2,0,1))

    projections = lsegs_tiled * dn  # points projected onto line segs
    distances = np.linalg.norm(point_segs-projections, axis=2)

    return distances.min(axis=1) # min distance from each point to all lines