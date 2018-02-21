import numpy as np

# def num(points):
#     expd=np.expand_dims(points,2) #need another dimension...
#     tiled=np.tile(expd, points.shape[0]) #...to tile up the vectors
#     trans=np.transpose(points) #Also need to transpose the points matrix to fit well with broadcasting
#     diff=trans-tiled           #doing the difference, exploiting Numpy broadcasting capabilities
#     num=np.sum(np.square(diff), axis=1) #an then obtain the squared norm of the difference
#     return num

# takes an NxMx2 array of point-line segments (point - start)
# and an Mx2 array of lines (end - start), returns NxM dot products
def elementwise_dot(points,lines):
    s = points.shape
    r = np.empty(s[0],s[1])
    for i in range(s0):
        for j in range(s1):
            r[i,j] = np.dot(points[i,j,:], lines[j,:])
    return r

# takes vectors of points and medial axis line segments, returns min distances
# using numpy's matrix operations for optimal performance
def dist_to_lines(points,lines):
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
    line_segs_diag = numpy.diag(line_segs)
    dots = ( np.matmul(point_segs(:,:,0),line_segs_diag(:,:,0)) + 
                 np.matmul(point_segs(:,:,1),line_segs_diag(:,:,1)) ) 
    #dots = elementwise_dot(point_segs, line_segs)
    
    norms = np.linalg.norm(line_segs, axis=1)
    dn = np.tile(np.expand_dims(dots/norms,2),2) # tile to match point_segs

    line_segs_exp = np.expand_dims(line_segs,2) # tile to match point_segs
    line_segs_tiled = np.tile(line_segs_exp, points.shape[0]).transpose((2,0,1))

    projections = line_segs_tiled * dn  # points projected onto line segs
    distances = np.linalg.norm(point_segs-projections, axis=2)

    return distances.min(axis=1).tolist()
