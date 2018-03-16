#!/usr/bin/python
 
import cv2
import numpy as np
import scipy.io as sio
import distances as dist
#from timeit import default_timer as timer



def printLines(lines) :
    for l in lines :
        cv2.line(img, tuple(l[0])[::-1], tuple(l[1])[::-1], (0,0,255))
    cv2.namedWindow('image', cv2.WINDOW_NORMAL)
    cv2.imshow('image', img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

def draw_point(img, p, color ) :
	cv2.circle( img, p, 1, color)#, cv2.FILLED, cv2.LINE_AA, 0 )

def printPoints(points) :
    for p in ima_points :
        cv2.circle( img, p, 1, (0,0,255) ) 
    cv2.namedWindow('image', cv2.WINDOW_NORMAL)
    cv2.imshow('image', img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
    
    
# Returns a set of 3x3 thinning elements which will thin skeleton to one-pixel
# width while preserving connectivity
def getThinningElements():
    x = 255
    struct_el = np.empty((8,3,3), dtype=np.uint8)
    
    struct_el[0] = np.array([[0, 0, 0],
                             [x, 1, x],
                             [1, 1, 1]])

    struct_el[2] = np.array([[0, x, 1],
                             [0, 1, 1],
                             [0, x, 1]])

    struct_el[4] = np.array([[1, 1, 1],
                             [x, 1, x],
                             [0, 0, 0]])

    struct_el[6] = np.array([[1, x, 0],
                             [1, 1, 0],
                             [1, x, 0]])

    struct_el[1] = np.array([[x, 0, 0],
                             [1, 1, 0],
                             [x, 1, x]])

    struct_el[3] = np.array([[0, 0, x],
                             [0, 1, 1],
                             [x, 1, x]])

    struct_el[5] = np.array([[x, 1, x],
                             [0, 1, 1],
                             [0, 0, x]])

    struct_el[7] = np.array([[x, 1, x],
                             [1, 1, 0],
                             [x, 0, 0]])
    return struct_el


# Returns a set of 3x3 thinning elements which will thin skeleton to one-pixel
# width while preserving connectivity
def getThinningElements2():
    x = 255
    struct_el = np.empty((16,3,3), dtype=np.uint8)
    
    struct_el[0] = np.array([[0, 0, 0],
                             [x, 1, x],
                             [1, 1, 1]])

    struct_el[4] = np.array([[0, x, 1],
                             [0, 1, 1],
                             [0, x, 1]])

    struct_el[8] = np.array([[1, 1, 1],
                             [x, 1, x],
                             [0, 0, 0]])

    struct_el[12] = np.array([[1, x, 0],
                              [1, 1, 0],
                              [1, x, 0]])

    struct_el[1] = np.array([[x, 0, 0],
                             [1, 1, 0],
                             [x, 1, x]])

    struct_el[5] = np.array([[0, 0, x],
                             [0, 1, 1],
                             [x, 1, x]])

    struct_el[9] = np.array([[x, 1, x],
                             [0, 1, 1],
                             [0, 0, x]])

    struct_el[13] = np.array([[x, 1, x],
                              [1, 1, 0],
                              [x, 0, 0]])
    
    struct_el[2] = np.array([[0, x, 1],
                             [0, 1, 1],
                             [0, x, x]])

    struct_el[6] = np.array([[1, 1, x],
                             [x, 1, x],
                             [0, 0, 0]])

    struct_el[10] = np.array([[x, x, 0],
                              [1, 1, 0],
                              [1, x, 0]])

    struct_el[14] = np.array([[0, 0, 0],
                              [x, 1, x],
                              [x, 1, 1]])
    
    struct_el[3] = np.array([[0, x, x],
                             [0, 1, 1],
                             [0, x, 1]])

    struct_el[7] = np.array([[x, 1, 1],
                             [x, 1, x],
                             [0, 0, 0]])

    struct_el[11] = np.array([[1, x, 0],
                              [1, 1, 0],
                              [x, x, 0]])

    struct_el[15] = np.array([[0, 0, 0],
                              [x, 1, x],
                              [1, 1, x]])
    
    return struct_el


# Calculates piecewise linear contour 
def approxContour (edge_points, contour_lines, contour_idx, thresh) :
    
    lines_to_split = []
    contour_dists = dist.pts2lines(edge_points,contour_lines)
    
    for i,l in enumerate(contour_lines) :
        contour_region = range(contour_idx[i],contour_idx[i+1])
        if max(contour_dists[contour_region,i]) > thresh :
            # calculate the index and contour location (i) of maximally distant point
            max_idx = np.argmax(contour_dists[contour_region,i])
            region_bound = max_idx+contour_idx[i]
            lines_to_split += [[i,region_bound]]
            
    if not lines_to_split:
        return contour_lines
    else :
        contour_lines_split = contour_lines
        for j in reversed(lines_to_split) : 
            # insert each maximally distant point into contour
            insert_from = np.array([contour_lines[j[0]][0], edge_points[j[1],:]])
            insert_to = np.array([edge_points[j[1],:], contour_lines[j[0]][1]])
            contour_lines_split[j[0]] = insert_from
            contour_lines_split = np.insert(contour_lines_split, j[0]+1, insert_to, axis=0)
            contour_idx = np.insert(contour_idx, j[0]+1, j[1])
        # apply recursively until all lines are within threshold value to original edge pts    
        return approxContour(edge_points, contour_lines_split, contour_idx, thresh)    


# Finds piecewise linear contour using helper function approxContour
def findApproxContour (img_bin, threshold) :
    
    # find points along shape edge (contour)
    edge_points = cv2.findContours(img_bin,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)[1][0]
    edge_points = np.flip(np.squeeze(edge_points, axis=1),axis=1)
    
    edge_min_y = np.where(edge_points[:,0]==np.min(edge_points[:,0]))
    min_y = edge_min_y[0][len(edge_min_y[0])/2]
    edge_max_y = np.where(edge_points[:,0]==np.max(edge_points[:,0]))
    max_y = edge_max_y[0][len(edge_max_y[0])/2]
    
    edge_min_x = np.where(edge_points[:,1]==np.min(edge_points[:,1]))
    min_x = edge_min_x[0][len(edge_min_x[0])/2]
    edge_max_x = np.where(edge_points[:,1]==np.max(edge_points[:,1]))
    max_x = edge_max_x[0][len(edge_max_x[0])/2]      
    
    c_idx = np.array([0,min_x,max_y,max_x,min_y,edge_points.shape[0]-1])
    
    c_lines = np.array([[edge_points[c_idx[0]], edge_points[c_idx[1]]],
                        [edge_points[c_idx[1]], edge_points[c_idx[2]]],
                        [edge_points[c_idx[2]], edge_points[c_idx[3]]],
                        [edge_points[c_idx[3]], edge_points[c_idx[4]]],
                        [edge_points[c_idx[4]], edge_points[c_idx[5]]]])
    
    # Recursively breaks contours into lines which are threshold-pixels close to the edge
    contour_lines = approxContour (edge_points, c_lines, c_idx, threshold)
    
    return contour_lines, edge_points


# Calclates the Integer Medial Axis for a given shape and feature-transform lines
def imaTransform(white_points,edge_points,ft_lines) :
     
    # Find all points in shape that are not boundary points
    interior_points = []
    for idx, point in enumerate(white_points):
        if np.all(np.any(edge_points!=point, axis=1)):
            interior_points.append([point.tolist(), idx])
    
    # Find all points which are adjacent (not diagonal) to points which have
    # a different closest contour line, these are the IMA point
    ima_points = []
    for ipoint in interior_points:
        closest_line = ft_lines[ipoint[1]]
        offsets = [[0,1],[0,-1],[1,0],[-1,0]]#[[0,1],[-1,0]]#
        inIMA = False
        for i in range(0,4):
            off_point = np.add(ipoint[0],offsets[i])
            off_line = ft_lines[np.where(np.all(white_points==off_point,axis=1))[0][0]]
            if closest_line != off_line:
                inIMA = True
        if inIMA:
            ima_points.append(ipoint[0])
            
    return ima_points


# Define the indices of the 8 neighboring pixels
def getNeighbors(pt) :
    offsets = [[0,0],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1],[1,0],[1,1]]
    return np.add(offsets,pt)


# Returns the number of zero-valued neighbors to point pt in binary img
def countNeighbors(pt, img_bin) :
    idxs = getNeighbors(pt)
    return np.sum(img_bin[idxs[:,0],idxs[:,1]]) - 1 # Subtract 1 for center pixel (pt) 


# Thins a skeleton by deleting any points which match the given struct element
def thin_element(points, img_bin, e) :
    hit_n_miss = np.zeros_like(img_ima_bin)
    pts2del = []
    # Scan each point with the struct element e and save any matching pixel
    for p in points:
        ps = getNeighbors(p)
        # Convert struct_el to P0-P8 format used by getNeighbors()
        struct_ps = np.array([e[1,1],e[1,2],e[0,2],
                              e[0,1],e[0,0],e[1,0],
                              e[2,0],e[2,1],e[2,2]])
        isMatch = True        
        # Test if each point in skeleton matches struct element pattern
        for i in range(9):
            if struct_ps[i]==255: 
                continue
            elif img_bin[ps[i,0],ps[i,1]] != struct_ps[i]:
                isMatch = False
        if isMatch:
            hit_n_miss[p[0],p[1]] = 1
            pts2del += [p]
    # Delete any match points from binary image and points list      
    img_bin = img_bin - hit_n_miss  
    pts2del = np.array(pts2del)
    for p in pts2del :
        idx = [np.where(np.all(points==p,axis=1))[0][0]]
        points = np.delete(points,idx,axis=0)
        
    return points, img_bin

    
# Repeat the hit-and-miss morphological thinning using struct_el until thinning
# operation returns the same shape
def thinSkeleton(skel_points, img_bin, struct_el) :
    thinned = False
    while not thinned:
        img_bin_prev = img_bin.copy()
        for i in range(struct_el.shape[0]):
            skel_points, img_bin = thin_element(skel_points, img_bin, struct_el[i])
        thinned = np.all(img_bin == img_bin_prev)  
    return skel_points, img_bin


# Remove any branches which are shorter than (thresh * [dist to closest contour line]) 
# Use values of thresh close to 1 to adjust pruning, 1.4 is a good first choice
def removeBranches(ima_points, img_ima_bin, threshold) :
    # Find end points of branches of the medial axis
    end_points = []
    for ima_p in ima_points:
        if countNeighbors(ima_p , img_ima_bin) < 2 :
            end_points += [ima_p ]
    end_points = np.array(end_points)

    # Find branches & delete those which fail to satisfy distance condition
    for end_p in end_points:
        branch_points = []
        prev = end_p
        curr = end_p
        while countNeighbors(curr, img_ima_bin) < 3 :
            # Identify next point along branch
            n = getNeighbors(curr)
            isNeighbor = (img_ima_bin[n[:,0],n[:,1]]!=0)
            isNeighbor[0] = False # neighbor 0 is the current pixel
            notPrev = np.any(n!=prev,axis=1)
            isNext = np.all([isNeighbor, notPrev],axis=0)
            # Add current point to branch, update to next point
            branch_points += [curr]
            prev = curr
            curr = np.squeeze(n[isNext])
        branch_points = np.array(branch_points)
        # Remove branches which are shorter than ft*thresh
        branch_len = dist.points2point([branch_points[-1]],branch_points[0])
        branch_ft = np.min(ft_dists[np.where(np.all(white_points==branch_points[-1],axis=1))[0][0]])
        if branch_len < branch_ft*threshold:
            for bp in branch_points :
                img_ima_bin[bp[0],bp[1]] = 0
                idx = [np.where(np.all(ima_points==bp,axis=1))[0][0]]
                ima_points = np.delete(ima_points,idx,axis=0)

    return ima_points, img_ima_bin


# Thin skeleton to 1-pixel width and remove any extraneous branches
def pruneSkeleton(ima_points, img_ima_bin, threshold) :
    # Thin IMA skeleton to one-pixel width
    struct_el = getThinningElements() # 1st pass elements
    ima_points, img_ima_bin = thinSkeleton(ima_points, img_ima_bin, struct_el) 
    
    # Prune branches on the IMA caused by joints in the contour line segments
    ima_points, img_ima_bin = removeBranches(ima_points, img_ima_bin, threshold)
     
    # Thin again to remove any remaining branch points
    struct_el = getThinningElements2() # 2nd pass - more aggressive elements
    ima_points, img_ima_bin = thinSkeleton(ima_points, img_ima_bin, struct_el) 
    
    return ima_points, img_ima_bin



if __name__ == '__main__':

    img_names = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12"]#["blake_01","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"] 
    pixel_threshold = [2,2,11,18,10,16,21,21]#,19] #,10,19,5,1,9,15]
    s_threshold = [1.4,1.2,11,18,10,16,21,21,19]
    patient = "DF"
    shape_img = ["Simple/" + img_name + ".png" for img_name in img_names]#["Blake/" + img_name + ".png" for img_name in img_names]
    out_path = './'+patient+'/shape_analysis/'
    
    super_matrix = zip(img_names,pixel_threshold,shape_img)

    # Determine pixel offset for contour detection
    offset = 10

    for row in super_matrix:
        #if row[0]!='solo3': continue
        print "Starting", row[0]
        
        # Read in the image.
        img = cv2.imread(row[2],cv2.IMREAD_UNCHANGED)

        # Convert alpha transparency to black
        img[(img[:,:,3]==0),0:3] = 0
        
        # Keep a copy around
        img_orig = img.copy();
        
        # Image must be offset to find contours
        img = cv2.copyMakeBorder(img, top=offset, bottom=offset, left=offset, right=offset, 
                                 borderType=cv2.BORDER_CONSTANT, value=0);
        
        # Binarize image
        img_bin = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
        img_bin[(img_bin!=0)] = 1 # convert to logical array (white = 1)
        
        # Compute image centroid
        white_count = np.sum(img_bin)
        white_idx = np.where(img_bin==1) # find all indices of white pixels
        centroid_y = int(np.sum(white_idx[0])/white_count)
        centroid_x = int(np.sum(white_idx[1])/white_count)
        centroid = [centroid_y,centroid_x,0]
        
        # Get piecewise linear approximation of shape edge which is threshold-pixels close
        # Allows us to perform the Integer Medial Axis transform with few unwanted branches
        contour_lines, edge_points = findApproxContour(img_bin, threshold=2)
       
        # Get list of all white points (shape interior)
        white_points = np.transpose(np.array(white_idx))

        # Calculate the feature transform for interior points
        ft_dists = dist.pts2lines(white_points,contour_lines) # distance from each white point to contour
        ft_lines = np.argmin(ft_dists,axis=1) # contour line closest to each white point

        # Find the Integer Medial Axis
        ima_points = imaTransform(white_points, edge_points, ft_lines)
        ima_points = np.array(ima_points)
    
        # Create a binary img representation of IMA points
        img_ima_bin = np.zeros_like(img_bin)
        img_ima_bin[ima_points[:,0],ima_points[:,1]] = 1
    
        # Prune branches on the IMA caused by joints in the contour line segments
        ima_points, img_ima_bin = pruneSkeleton(ima_points, img_ima_bin, threshold=1.2)

        ima_points, img_ima_bin = pruneSkeleton(ima_points, img_ima_bin, threshold=0.2)

        # Format data to x-y coordinates (OpenCV uses y-x), remove offset
        ima_points = np.flip(ima_points, axis=1) - [offset,offset]
        edge_points = np.flip(edge_points, axis=1) - [offset,offset]
        centroid = [centroid[1]-offset, centroid[0]-offset]
        
        # Draw image with medial axis points
        for p in ima_points :
            cv2.circle( img_orig, tuple(p), 1, (255,0,0,255) ) 
        for e in edge_points :
            cv2.circle( img_orig, tuple(e), 1, (0,0,255,255) ) 
        cv2.circle( img_orig, tuple(centroid), 4, (0,255,0,255) )    

#       # DEBUG - show image
#        cv2.namedWindow('image', cv2.WINDOW_NORMAL)
#        cv2.imshow('image', img_orig)
#        cv2.waitKey(0)
#        cv2.destroyAllWindows()
#        break
            
#ima_pts = np.where(img_ima_bin==1)
#ima_pts = np.array(ima_pts)
#ima_pts = np.transpose(ima_pts)
#ima_pts = np.flip(ima_pts, axis=1) - [offset,offset]
#ima_points = ima_pts
#img = cv2.imread(row[2],cv2.IMREAD_UNCHANGED)
## Convert alpha transparency to black
#img[(img[:,:,3]==0),0:3] = 0
## Keep a copy around
#img_orig = img.copy();
#for p in ima_points :
#    cv2.circle( img_orig, tuple(p), 1, (255,0,0,255) ) 
#for e in edge_points :
#    cv2.circle( img_orig, tuple(e), 1, (0,0,255,255) ) 
#cv2.circle( img_orig, tuple(centroid), 4, (0,255,0,255) ) 
#
#cv2.namedWindow('image', cv2.WINDOW_NORMAL)
#cv2.imshow('image', img_orig)
#cv2.waitKey(0)
#cv2.destroyAllWindows()

        # Save a copy of the image
        cv2.imwrite(out_path + row[0] + '_ma.png',img_orig)
        
        #Save Medial Axis, Edge Points and Centroid to MAT file
        sio.savemat(out_path + row[0] + '_shape_analysis.mat', {'ma_points':ima_points, 'edge_points':edge_points, 'centroid':centroid});


