function canvas = composeImage(events,img_dir,canvas_width, canvas_height)
    image_files = getAllFiles(img_dir);
    canvas = [canvas_height, canvas_width];
    if ~isempty(image_files) % we do have images and we're not crazy.
        % look for the image that we want
        for i = 1:size(events,1)
            for j = 1:size(image_files,1)
                % correct image, overlay it
                if contains(image_files{j},events{i,1}{1})
                    % grab position and rotation data
                    [event,map,alpha] = imread(image_files{j},'png');
                    scalefactor = events{i,1}{5}/canvas_height;
                    
                    
                    % 4 is rotation
                    imrotate(event,events{i,1}{4});
                    imrotate(alpha,events{i,1}{4});
                    
                    
                    
                    % rotate image, place centred at position (percentage
                    % of screen, need to check how unity handles)
                end
            end
        end
    end
end