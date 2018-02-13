#!/usr/bin/python
 
import cv2
import numpy as np
import random
import pdb


def isWhite(point):
	y=point[0]
	x=point[1]
	
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
	shape_img = "images/blake_12.png"
	divvy = 19


	shape_anal = True
	path = "./results_trial_scaled_absolute/vPzW_results.txt"

 
	# Define window names
	win_delaunay = "Delaunay Triangulation"
	win_voronoi = "Voronoi Diagram"
 
	# Turn on animation while drawing triangles
	animate = False
	 
	# Define colors for drawing.
	delaunay_color = (255,255,255)
	points_color = (0, 0, 255)
 
	# Read in the image.
	img = cv2.imread(shape_img,cv2.IMREAD_COLOR);
	 
	# Keep a copy around
	img_orig = img.copy();
	 
	# Rectangle to be used with Subdiv2D
	size = img.shape
	rect = (0, 0, size[1], size[0])
	 
	# Create an instance of Subdiv2D
	subdiv = cv2.Subdiv2D(rect);

	# Create an array of points.
	point = [size[0]/2,1]
	points = []
	black = [size[0]/2,0]

	whites = [0, 0]
	white_count = 0 #################################################################################################### make this a list of shapes - white counts
	for x in range(0,size[0]):
		for y in range(0,size[1]):
			if isWhite([x, y]):
				whites[0] += x
				whites[1] += y
				white_count += 1
	print white_count

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
		if count % divvy == 0:
			temp_points += [tuple([point[1],point[0]])]

	points = temp_points

	# pdb.set_trace()


	# Insert points into subdiv
	myCnt = []
	for p in points :
		myCnt += [[p[0],p[1]]]
		subdiv.insert(p)

 
	# Draw delaunay triangles
	draw_delaunay( img, subdiv, (255, 255, 255) );
	# Draw points
	for p in points :
		draw_point(img, p, (0,0,255))
 
	# Allocate space for Voronoi Diagram
	img_voronoi = np.zeros(img.shape, dtype = img.dtype)
 
	# Draw Voronoi diagram
	ma_lines = draw_voronoi(img_voronoi,subdiv,myCnt)
    print(ma_lines)



	import csv
	import assist

	with open(path) as f:
	    reader = csv.reader(f, delimiter="\t")
	    d = list(reader)
	this_shape = []
	for s in d:
		if (s[0]+".png") == shape_img.split('/')[-1]:
			s_x = float(s[1])
			s_y = (float(s[2])-size[0])*-1
			this_shape += [[s_x, s_y, 0]]



	# pdb.set_trace()

	min_dists = []
	for touch in this_shape:
		min_dist = 9999
		for line in ma_lines:
			dist = assist.pnt2line(touch, np.append(line[0], 0), np.append(line[1], 0))[0]

			if dist < min_dist:
				min_dist = dist


		min_dists += [min_dist]

	# text1 = 'Center (X, Y): '+ str(round(center_x, 2)) + ', ' + str(round(center_y, 2))
	text2 = 'MSE Axis: ' + str(round(np.average(np.power(min_dists, 2)), 2))
	# text3 = 'MSE Centre: ' + str(round(np.average(np.power(center_dists, 2)), 2))
	fact = 'Branch Factor: ' + str(divvy)

	# cv2.putText(img_voronoi, text1, (10,20), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 0, 1)
	# cv2.putText(img_voronoi, text2, (10,50), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 0, 1)
	# cv2.putText(img_voronoi, text3, (10,80), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 0, 1)
	# cv2.putText(img_voronoi, fact, (10,110), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 0, 1)

	# Show results
	cv2.imshow(win_delaunay,img)
	cv2.imshow(win_voronoi,img_voronoi)
	cv2.waitKey(0)


	# print 'Shape\tCentreX\tCenterY\tMSE-axis\tMSE-center\tBranching Factor'
	print shape_img[:-4] + '\t'+ str(round(np.average(np.power(min_dists, 2)), 2))+'\t' + str(divvy)



# 	[array([  68.85739899,  393.21841431], dtype=float32), array([  87.0056839 ,  396.32955933], dtype=float32)]
# number of lines =  320
# length of lines =  5336.76109763
# shape10	626.86	391.64	2348.54	24393.31	18