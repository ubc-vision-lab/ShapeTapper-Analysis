% import cv2
% import numpy as np
% import scipy.io as sio
% import vectorization as vect
% 
% import assist

global nextCirc
nextCirc = [ [-1,-1]; [0,1]; [1,-1]; [1,0];...
             [1,1];   [0,1]; [-1,1]; [-1,0]; ];

global img
global img_size

% Define colors for drawing.
delaunay_color = [255,255,255];
points_color = [0, 0, 255];
    
img_names = ['blake_01';'blake_04';'blake_06';'blake_07';'blake_08';'blake_09';'blake_10';'blake_11';'blake_12']; %,"solo10","solo11","solo3","solo6","solo7","solo9"]
divvy = [18,22,11,18,10,16,21,21,19]; %,10,19,5,1,9,15]
patient = 'DF';
%shape_img = strcat('Blake/',strcat(img_names,'.png'));
%generated_files = strcat(['./' patient '/generated_results/1000/'],strcat(img_names,['_Patient_' patient '_generated_results.mat']));
%observed_files = strcat(['./' patient '/aggregated_observations/'],strcat(img_names,['_Patient_' patient '_aggregated_observations.mat']));
out_path = [patient '/test_mat/'];
%print observed_files
                
%super_matrix = zip(img_names,divvy,shape_img,generated_files,observed_files);
                
aggregate_ma = [];
aggregate_c = [];
                
%global_start = timer()
for f=1:length(img_names)
    
    %start_row = timer()
    
    % Read in the image.
    [img,map] = imread(['Blake/' img_names(f,:) '.png']);
    img = ind2rgb(img, map);
    img_orig = img; % Keep a copy around

    % Rectangle to be used with Subdiv2D
    img_size = size(img);
    rect = [0, 0, img_size(2), img_size(1)];

    % Create an instance of Subdiv2D
    %subdiv = cv2.Subdiv2D(rect);

    % Create an array of points.
    point = [round(img_size(1)/3),1];
    points = [];
    black = [round(img_size(1)/3),0];

    whites = [0,0];
    white_count = 0;     % make this a list of shapes - white counts
    for y=1:img_size(1)     % img_size(1) is number of rows
        for x=1:img_size(2) % img_size(2) is number of columns
            if isWhite([y,x])
                whites = whites + [y,x];
                white_count = white_count + 1;
            end
        end
    end
    centroid = [whites/white_count, 0];

    while 1 % go until you get to a white pixel
        if isWhite(point)
            points = [points; point];
            break
        else
            black = point;
            point(2) = point(2) + 1;
        end
    end

    while 1 % find all
        [black, nextPoint] = nextPixel(black, point, points);
        if all(nextPoint~=[-1,-1]) && ~isWhite(black)
            point = nextPoint;
            points = [points; point];
        else
            break
        end
    end
    
    count = 0;
    temp_points = [];
    for p=1:divvy(f):length(points)
        temp_points = [temp_points; points(p,:)];
    end
    points = temp_points;
    
    % Insert points into subdiv
    myCnt = [];
    for p=1:length(points)
        myCnt = [myCnt; points(p,:)];
        %subdiv.insert(points(p,:));
    end
    
    % Draw delaunay triangles
    %draw_delaunay( img, subdiv, delaunay_color );
    % Draw points
%     for p=1:length(points)
%         draw_point(img, points(p,:), points_color);
%     end

    % Allocate space for Voronoi Diagram
    img_voronoi = zeros(size(img));

    % Draw Voronoi diagram
    ma_lines = draw_voronoi(img_voronoi,subdiv,myCnt);
    % sio.savemat('./Patient_'+patient+'/ma_c_values/' + row[0] + '_ma.mat',{'ma_lines':ma_lines});
    % sio.savemat('./Patient_'+patient+'/ma_c_values/' + row[0] + '_c.mat',{'centroid':centroid});

    %print "Calculating Observed"
    % let's start computing distances!
    % observed first. Faster.
    %observed_mat = sio.loadmat(row(4));
    %observed = observed_mat['img_dataset']

    %start_time = timer()
    observed_min_dists = vect.dist_pts2lines(observed,ma_lines);
    observed_centroid_dists = vect.dist_pts2cent(observed,centroid);
    %print "Observed: ",(timer() - start_time)*1000, 'ms'

    observed_medaxis_var = sum(power(observed_min_dists, 2))/(size(observed_min_dists)-1);
    observed_medaxis_amd = mean(observed_min_dists);
    observed_medaxis_data = [observed_medaxis_var, observed_medaxis_amd];
%     sio.savemat(out_path+row(0)+'_medaxis_new.mat',{'observed_medaxis_data':observed_medaxis_data});

    observed_centroid_var = sum(power(observed_centroid_dists, 2))/(size(observed_centroid_dists)-1);
    observed_centroid_amd = mean(observed_centroid_dists);
    observed_centroid_data = [observed_centroid_var, observed_centroid_amd];
%     sio.savemat(out_path+row(0)+'_centroid_new.mat',{'observed_centroid_data':observed_centroid_data});

    %print "Calculating Generated"
    %generated_mat = sio.loadmat(row(3));
    %generated = generated_mat['img_datasets'];
    generated_centroid_data = [];
    generated_medaxis_data = [];

    generated_min_dists = zeros(generated.shape(0),generated.shape(1));
    generated_centroid_dists = zeros(generated.shape(0),generated.shape(1));
    %start_time = timer()
    for i=1:length(generated)%, point_set in enumerate(generated)
        generated_min_dists(i) = vect.dist_pts2lines(point_set,ma_lines);
        generated_centroid_dists(i) = vect.dist_pts2cent(point_set,centroid);
        %print "Generated: ",(timer() - start_time), 's'
    end
    % medial axis variables
    generated_medaxis_var = sum(np.power(generated_min_dists,2)) / (generated_min_dists.shape(1) - 1);
    generated_medaxis_amd = mean(generated_min_dists);%, axis=1)
    generated_medaxis_data = [generated_medaxis_var, generated_medaxis_amd];

    generated_centroid_var = sum(np.power(generated_centroid_dists,2)) / (generated_centroid_dists.shape(1) - 1);
    generated_centroid_amd = mean(generated_centroid_dists);
    generated_centroid_data = [generated_centroid_var, generated_centroid_amd];
end
    % print "saving file"
%     generated_medaxis_data = [list(datum) for datum in generated_medaxis_data]
%     generated_centroid_data = [list(datum) for datum in generated_centroid_data]
% 
%     sio.savemat(out_path+row[0]+'_generated_ma_new.mat',{'generated_medaxis_data':generated_medaxis_data})
%     sio.savemat(out_path+row[0]+'_generated_c_new.mat',{'generated_centroid_data':generated_centroid_data})

                
function w = isWhite(point)
    global img
    global img_size
    w = 0;
    if any(point > img_size(1:2)) || any(point <= 0)
        return
    elseif all(squeeze(img(point(1),point(2),:)))
        w = 1;
    end
end                                                                                
                                                                                
% finds the first white point that is not in points
function [black,newPoint] = nextPixel(black, point, points)
    global nextCirc
    for w=1:8
        pt_offset = point + nextCirc(w,:);
        % if offset point is white and not already in points
        if isWhite(pt_offset) && ~ismember(2,sum(points==pt_offset,2))
            newPoint = pt_offset;    
            return
        else
            black = pt_offset;
        end
    end
    newPoint = [-1,-1];
    return
end 

% Check if a point is inside a rectangle
function contains = rect_contains(rect, point)
    contains = 1;
    if any(point<rect(1:2)) || any(point>rect(3:4))
        contains = 0;
    end
end

% Draw a point
function draw_point(img, p, color )
    cv2.circle( img, p, 2, color, cv2.FILLED, cv2.LINE_AA, 0 );
end
 
%  
% % Draw delaunay triangles
% def draw_delaunay(img, subdiv, delaunay_color ) :
%  
%     triangleList = subdiv.getTriangleList();
%     img_size = img.shape
%     r = (0, 0, img_size[1], img_size[0])
%  
%     for t in triangleList :
%          
%         pt1 = (t[0], t[1])
%         pt2 = (t[2], t[3])
%         pt3 = (t[4], t[5])
%          
%         if rect_contains(r, pt1) and rect_contains(r, pt2) and rect_contains(r, pt3) :
%          
%             cv2.line(img, pt1, pt2, delaunay_color, 1, cv2.LINE_AA, 0)
%             cv2.line(img, pt2, pt3, delaunay_color, 1, cv2.LINE_AA, 0)
%             cv2.line(img, pt3, pt1, delaunay_color, 1, cv2.LINE_AA, 0)
%  
%  
% % Draw voronoi diagram
% def draw_voronoi(img, subdiv, myCnt) :
%     img_size = img.shape
%     rect = np.array([[0,0], [0,img_size[0]], [img_size[1],img_size[0]], [img_size[1],0]])
%     cv2.fillConvexPoly(img, rect, (255,255,255), cv2.LINE_AA, 0)
% 
%     [facets, centers] = subdiv.getVoronoiFacetList([])
%     lines = []
%     
%     for i in xrange(0,len(facets)) :
%         for f in range(0,len(facets[i])) :
%             if cv2.pointPolygonTest(np.array(myCnt), tuple(facets[i][f]), False) > 0 && cv2.pointPolygonTest(np.array(myCnt), tuple(facets[i][(f+1)%len(facets[i])]), False) > 0 :
%                 lines.append([facets[i][f], facets[i][(f+1);%len(facets[i])]])
%             end
%         cv2.circle(img, (centers[i][0], centers[i][1]), 3, (0, 0, 0), cv2.FILLED, cv2.LINE_AA, 0);
%         end
%     end
%         
%     for l in lines:
%         cv2.line(img, tuple(l[0]), tuple(l[1]), (0,0,0), 1, cv2.LINE_AA)
%         % print l[0][0], l[0][1], l[1][0], l[1][1]
%     return lines
