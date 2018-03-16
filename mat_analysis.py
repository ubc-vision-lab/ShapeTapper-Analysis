#!/usr/bin/python
 
import cv2
import numpy as np
import scipy.io as sio
import distances as dist
#from timeit import default_timer as timer


def getVarMean(data) :
    var = np.sum(np.power(data, 2))/(np.size(data)-1)
    amd = np.mean(data)
    return [var, amd]

def getVarMean2D(data) :
    var = np.sum(np.power(data,2),axis=1) / (data.shape[1] - 1)
    amd = np.mean(data, axis=1)
    return [var, amd]


if __name__ == '__main__':

    img_names = ["blake_01","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"] #,"solo10","solo11","solo3","solo6","solo7","solo9"]
    patient = "DF"
    shape_img = ["Blake/" + img_name + ".png" for img_name in img_names]
    observed_files = ["./"+patient+"/aggregated_observations/" + img_name + "_Patient_"+patient+"_aggregated_observations.mat" for img_name in img_names]
    generated_files = ["./"+patient+"/generated_results/1000/" + img_name + "_Patient_"+patient+"_generated_results.mat" for img_name in img_names]
    shape_files = ["./"+patient+"/shape_analysis/" + img_name + "_shape_analysis.mat" for img_name in img_names]
    out_path = "./"+patient+"/shape_analysis/"
    
    super_matrix = zip(img_names,shape_img,shape_files,observed_files,generated_files)

    for row in super_matrix:
        #if row[0]!='blake_01': continue
        print "Starting", row[0]
        #start_row = timer()
        #print row      
        
        # Read in the image.
        img = cv2.imread(row[1],cv2.IMREAD_UNCHANGED)
        
        shape_analysis = sio.loadmat(row[2])
        ma_points = shape_analysis['ma_points']
        edge_points = shape_analysis['edge_points']
        centroid = shape_analysis['centroid']

        #print "Calculating Observed"
        observed_mat = sio.loadmat(row[3])
        observed = observed_mat['img_dataset']
        observed[:,1] = img.shape[0]-observed[:,1] # opencv coordinates use inverted y-axis

        #start_time = timer()
        observed_ma_min_dists = np.min(dist.points2points(observed,ma_points),axis=1)
        observed_edge_min_dists = np.min(dist.points2points(observed,edge_points),axis=1)
        observed_centroid_dists = dist.points2point(observed,centroid)
        #print "Observed: ",(timer() - start_time)*1000, 'ms'
        
        observed_medaxis_data = getVarMean(observed_ma_min_dists)
        observed_edge_data = getVarMean(observed_edge_min_dists)
        observed_centroid_data = getVarMean(observed_centroid_dists)
        sio.savemat(out_path+row[0]+'_observed_analysis.mat',
                    {'observed_medaxis_data':observed_medaxis_data, 
                     'observed_edge_data':observed_edge_data, 
                     'observed_centroid_data':observed_centroid_data})

        #print "Calculating Generated"
        generated_mat = sio.loadmat(row[4])
        generated = generated_mat['img_datasets']

        generated_ma_min_dists = np.empty((generated.shape[0],generated.shape[1]))
        generated_edge_min_dists = np.empty((generated.shape[0],generated.shape[1]))
        generated_centroid_dists = np.empty((generated.shape[0],generated.shape[1]))
        #start_time = timer()
        for i, point_set in enumerate(generated) :
            generated_ma_min_dists[i] = np.min(dist.points2points(point_set,ma_points),axis=1)
            generated_edge_min_dists[i] = np.min(dist.points2points(point_set,edge_points),axis=1)
            generated_centroid_dists[i] = dist.points2point(point_set,centroid)
        #print "Generated: ", (timer() - start_time), 's'
#        
#        # print "saving file"
        generated_medaxis_data = getVarMean2D(generated_ma_min_dists)
        generated_edge_data = getVarMean2D(generated_edge_min_dists)
        generated_centroid_data = getVarMean2D(generated_centroid_dists)
        generated_medaxis_data = [list(datum) for datum in generated_medaxis_data]
        generated_edge_data = [list(datum) for datum in generated_edge_data]
        generated_centroid_data = [list(datum) for datum in generated_centroid_data]
        sio.savemat(out_path+row[0]+'_generated_analysis.mat',
                    {'generated_medaxis_data':generated_medaxis_data, 
                     'generated_edge_data':generated_edge_data,
                     'generated_centroid_data':generated_centroid_data})
        
        # Draw image with medial axis points
        img[(img[:,:,3]==0),0:3] = 0 # Convert alpha transparency to black
        for p in ma_points :
            cv2.circle( img, tuple(p), 1, (255,0,0,255), thickness=1 ) 
        for o in observed :
            cv2.circle( img, tuple([int(o[0]),int(o[1])]), 1, (0,0,255,255), thickness=-1 ) 
        cv2.circle( img, tuple(np.squeeze(centroid).tolist()), 4, (0,255,0,255), thickness=2)
        
#        cv2.namedWindow('image', cv2.WINDOW_NORMAL)
#        cv2.imshow('image', img)
#        cv2.waitKey(0)
#        cv2.destroyAllWindows()
#        break
    
        cv2.imwrite(out_path + row[0] + '_touch_data.png',img)

