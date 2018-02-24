#!/usr/bin/python
 
import cv2
import numpy as np
import random
import pdb
import scipy.io as sio
import csv
import assist
import vectorization

def isWhite(point):
	y=point[0]
	x=point[1]
	if y>=img.shape[0] or x >=img.shape[1]: return False
	if y<0 or x <0: return False
	return img[y,x][0] == 255 and img[y,x][1] == 255 and img[y,x][2] == 255

def nextCirc(y,x):

	if y == -1 and x == -1:
		y = 0
		x = -1
	elif y == 0 and x == -1:
		y = 1
		x = -1
	elif y == 1 and x == -1:
		y = 1
		x = 0
	elif y == 1 and x == 0:
		y = 1
		x = 1
	elif y == 1 and x == 1:
		y = 0
		x = 1
	elif y == 0 and x == 1:
		y = -1
		x = 1
	elif y == -1 and x == 1:
		y = -1
		x = 0
	else:
		y = -1
		x = -1


	return [y,x]


# finds the first white point that is not in points
def nextpixel(black, point, points):
	y=point[0]
	x=point[1]

	for w in range(0,8):
		offset = nextCirc(black[0]-point[0], black[1]-point[1])
		if isWhite([y+offset[0],x+offset[1]]) and [y+offset[0],x+offset[1]] not in points:
			newPoint = [y+offset[0],x+offset[1]]	
			return black, newPoint
		else:
			black = [y+offset[0],x+offset[1]]
	print False
	return black, False
 

# Check if a point is inside a rectangle
def rect_contains(rect, point) :
	if point[0] < rect[0] :
		return False
	elif point[1] < rect[1] :
		return False
	elif point[0] > rect[2] :
		return False
	elif point[1] > rect[3] :
		return False
	return True
 

# Draw a point
def draw_point(img, p, color ) :
	cv2.circle( img, p, 2, color, cv2.FILLED, cv2.LINE_AA, 0 )
 
 
# Draw delaunay triangles
def draw_delaunay(img, subdiv, delaunay_color ) :
 
	triangleList = subdiv.getTriangleList();
	size = img.shape
	r = (0, 0, size[1], size[0])
 
	for t in triangleList :
		 
		pt1 = (t[0], t[1])
		pt2 = (t[2], t[3])
		pt3 = (t[4], t[5])
		 
		if rect_contains(r, pt1) and rect_contains(r, pt2) and rect_contains(r, pt3) :
		 
			cv2.line(img, pt1, pt2, delaunay_color, 1, cv2.LINE_AA, 0)
			cv2.line(img, pt2, pt3, delaunay_color, 1, cv2.LINE_AA, 0)
			cv2.line(img, pt3, pt1, delaunay_color, 1, cv2.LINE_AA, 0)
 
 
# Draw voronoi diagram
def draw_voronoi(img, subdiv, myCnt) :
	size = img.shape
	rect = np.array([[0,0], [0,size[0]], [size[1],size[0]], [size[1],0]])
	cv2.fillConvexPoly(img, rect, (255,255,255), cv2.LINE_AA, 0)

	( facets, centers) = subdiv.getVoronoiFacetList([])
	lines = []
	

	
	for i in xrange(0,len(facets)) :
		for f in range(0,len(facets[i])) :
			if cv2.pointPolygonTest(np.array(myCnt), tuple(facets[i][f]), False) > 0 and cv2.pointPolygonTest(np.array(myCnt), tuple(facets[i][(f+1)%len(facets[i])]), False) > 0 :
				lines.append([facets[i][f], facets[i][(f+1)%len(facets[i])]])




		cv2.circle(img, (centers[i][0], centers[i][1]), 3, (0, 0, 0), cv2.FILLED, cv2.LINE_AA, 0)
	
	for l in lines:
		cv2.line(img, tuple(l[0]), tuple(l[1]), (0,0,0), 1, cv2.LINE_AA)
		# print l[0][0], l[0][1], l[1][0], l[1][1]
	return lines




if __name__ == '__main__':

	img_names = ["blake_01","blake_04"]#,"blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"] #,"solo10","solo11","solo3","solo6","solo7","solo9"]
	# img_names = ["solo11"]
	divvy = [18,22]#,11,18,10,16,21,21,19] #,10,19,5,1,9,15]
	# divvy = [19]
	patient = "DF"
	shape_img = ["Blake/" + img_name + ".png" for img_name in img_names]
	generated_files = ["./"+patient+"/generated_results/1000/" + img_name + "_Patient_"+patient+"_generated_results.mat" for img_name in img_names]
	observed_files = ["./"+patient+"/aggregated_observations/" + img_name + "_Patient_"+patient+"_aggregated_observations.mat" for img_name in img_names]
	
	print observed_files

 	super_matrix = zip(img_names,divvy,shape_img,generated_files,observed_files)

 	aggregate_ma = [];
 	aggregate_c = [];
 	for row in super_matrix:
 		print row
 		# Define colors for drawing.
		delaunay_color = (255,255,255)
		points_color = (0, 0, 255)
	 
		# Read in the image.
		img = cv2.imread(row[2],cv2.IMREAD_UNCHANGED);
		 
		# Keep a copy around
		img_orig = img.copy();

		# Convert alpha transparency to black
		img[(img[:,:,3]==0),0:3] = 0

		# Rectangle to be used with Subdiv2D
		size = img.shape
		rect = (0, 0, size[1], size[0])
		 
		# Create an instance of Subdiv2D
		subdiv = cv2.Subdiv2D(rect);

		# Create an array of points.
		point = [size[0]/3,1]
		points = []
		black = [size[0]/3,0]

		whites = [0, 0]
		white_count = 0 # make this a list of shapes - white counts
		for y in range(0,size[0]): # size[0] is number of rows
			for x in range(0,size[1]): # size[1] is number of columns
				# consequently, size[2] is RGB
				if isWhite([y, x]):
					whites[0] += y
					whites[1] += x
					white_count += 1

		centroid_y = whites[0]/white_count
		centroid_x = whites[1]/white_count
		centroid = [centroid_y, centroid_x,0]

		while True: # go until you get to a white pixel
			if isWhite(point):
				points += [point]
				break
			else:
				black = point
				point[1] += 1
		
		while True: # find all 
			black, nextPoint =  nextpixel(black, point, points)

			if nextPoint and not isWhite(black):
				point = nextPoint
				points += [point]
			else:
				break

		count = 0
		temp_points = []
		for point in points:
			count += 1
			if count % row[1] == 0:
				temp_points += [tuple([point[1],point[0]])]

		points = temp_points

		# Insert points into subdiv
		myCnt = []
		for p in points :
			myCnt += [[p[0],p[1]]]
			subdiv.insert(p)

	 	print "Drawing Delaunay"
		# Draw delaunay triangles
		draw_delaunay( img, subdiv, (255, 255, 255) );
		# Draw points
		for p in points :
			draw_point(img, p, (0,0,255))
	 
		# Allocate space for Voronoi Diagram
		img_voronoi = np.zeros(img.shape, dtype = img.dtype)
	 	
	 	print "Drawing Voronoi"
		# Draw Voronoi diagram
		ma_lines = draw_voronoi(img_voronoi,subdiv,myCnt)
		print np.shape(ma_lines)
		ma_lines = np.array(ma_lines)

		# sio.savemat('./Patient_'+patient+'/ma_c_values/' + row[0] + '_ma.mat',{'ma_lines':ma_lines});
		# sio.savemat('./Patient_'+patient+'/ma_c_values/' + row[0] + '_c.mat',{'centroid':centroid});

		print "Calculating Observed"
		# let's start computing distances!
		# observed first. Faster.
		observed_mat = sio.loadmat(row[4])
		observed = observed_mat['img_dataset']
		observed_centroid_data = []
		observed_medaxis_data = []
		print np.shape(observed)

		observed_min_dists = vectorization.dist_pts2lines(observed,ma_lines)
		print np.shape(observed_min_dists)

		# # distance utility function is three dimensional, so fill in the 'z' axis with zeros
		# filler = np.transpose([np.zeros(np.shape(observed)[0],dtype=int)]) # works as advertised
		# filler works for both as the number of elements between the datasets is the same, but we can always make it
		# recalculated if needed
		# observed_filled = np.concatenate((observed,filler),axis=1)
		# # observed_min_dists = []
		# # observed_centroid_dists = []
		# # for observation in observed_filled:
		# # 	min_dist = 999999
		# # 	observed_centroid_dists += [assist.distance(observation,centroid)]
		# # 	for line in ma_lines:
		# # 		dist = assist.pnt2line(observation, np.append(line[0], 0), np.append(line[1], 0))[0]

		# # 		if dist < min_dist:
		# # 			min_dist = dist
		# # 	observed_min_dists += [min_dist]

		# # observed_medaxis_var = np.sum(np.power(observed_min_dists, 2))/(np.size(observed_min_dists)-1)
		# # observed_medaxis_amd = np.mean(observed_min_dists)
		# # observed_medaxis_data = [observed_medaxis_var, observed_medaxis_amd]
		# # sio.savemat(row[0]+'medaxis.mat',{'observed_medaxis_data':observed_medaxis_data})

		# # observed_centroid_var = np.sum(np.power(observed_centroid_dists, 2))/(np.size(observed_centroid_dists)-1)
		# # observed_centroid_amd = np.mean(observed_centroid_dists)
		# # observed_centroid_data = [observed_centroid_var, observed_centroid_amd]
		# # sio.savemat(row[0]+'centroid.mat',{'observed_centroid_data':observed_centroid_data})


		print "Calculating Generated"
		generated_mat = sio.loadmat(row[3])
		generated = generated_mat['img_datasets']
		generated_centroid_data = []
		generated_medaxis_data = []
		print np.shape(generated)

		generated_min_dists = vectorization.dist_pts2lines(generated,ma_lines)
		print np.shape(generated_min_dists)

		# count = 0
		# # this loop is for the generated data
		# filler = np.transpose([np.zeros(np.shape(generated)[1],dtype=int)]) # works as advertised
		# print np.shape(generated)
		# for result in generated:
		# 	result_filled = np.concatenate((result,filler),axis=1)
		# 	count += 1
		# 	generated_min_dists = []
		# 	generated_centroid_dists = []
		# 	if count % 10000 == 0:
		# 		print "getting distributions"
		# 	for touch in result_filled:
		# 		min_dist = 9999
		# 		generated_centroid_dists += [assist.distance(touch,centroid)]
		# 		for line in ma_lines:
		# 			dist = assist.pnt2line(touch, np.append(line[0], 0), np.append(line[1], 0))[0]

		# 			if dist < min_dist:
		# 				min_dist = dist


		# 		generated_min_dists += [min_dist]

		# 	# medial axis variables
		# 	generated_medaxis_var = np.sum(np.power(generated_min_dists, 2))/(np.size(generated_min_dists)-1)
		# 	# medaxis_mse = np.mean(np.power(generated_min_dists, 2))
		# 	# medaxis_std = np.sqrt(medaxis_var)
		# 	generated_medaxis_amd = np.mean(generated_min_dists)
		# 	generated_medaxis_data.append((generated_medaxis_var, generated_medaxis_amd))

		# 	generated_centroid_var = np.sum(np.power(generated_centroid_dists, 2))/(np.size(generated_centroid_dists)-1)
		# 	generated_centroid_amd = np.mean(generated_centroid_dists)
		# 	generated_centroid_data.append((generated_centroid_var, generated_centroid_amd))

		# print "saving file"
		# generated_medaxis_data = [list(datum) for datum in generated_medaxis_data]

		# generated_centroid_data = [list(datum) for datum in generated_centroid_data]
		# sio.savemat('./Patient_'+patient+'/'+ row[0]+'_ma.mat',{'generated_medaxis_data':generated_medaxis_data})
		# sio.savemat('./Patient_'+patient+'/'+ row[0]+'_c.mat',{'generated_centroid_data':generated_centroid_data})
		# # print generated_medaxis_data
		# # print np.shape(generated_medaxis_data)
		# # print generated_centroid_data
		# # print np.shape(generated_centroid_data)