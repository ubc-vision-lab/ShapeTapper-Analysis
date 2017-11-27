function [image_overlay, alpha_mask] = composeImage(events,img_dir,screen_width, screen_height, dpi)
    x_padding = 300;
    y_padding = 300;
    image_files = getAllFiles(img_dir);
    image_overlay = zeros(screen_height+y_padding, screen_width+x_padding, 3, 'uint8');
    alpha_mask = zeros(screen_height+y_padding, screen_width+x_padding, 'uint8');
    if ~isempty(image_files) % we do have images and we're not crazy.
        % look for the image that we want
        for i = 1:size(events,1)
            for j = 1:size(image_files,1)
                % correct image, overlay it
                if contains(image_files{j},events{i,1}{1})
                    % grab position and rotation data
                    [event,~,alpha] = imread(image_files{j},'BackgroundColor','none');
                    diag = norm([size(event,1),size(event,2)]);
                    diag_max = (screen_height * events{i,1}{5}) / 100;
                    scale = double(diag_max) / diag;
                    
                    % scale and rotate event, alpha channel
                    resized_event = imresize(event,scale);
                    resized_rotated_event = imrotate(resized_event,events{i}{4});
                    resized_alpha = imresize(alpha,scale);
                    resized_rotated_alpha = imrotate(resized_alpha,events{i}{4});
                    
                    % put the image and its alpha channel onto the white
                    % background
                    [resized_rotated_y, resized_rotated_x,~] = size(resized_rotated_event);
                    
                    world_to_pixels = screen_height/5/2; % emulating orthographic size (set to 5 in unity);
                    border = 0;%dpi/world_to_pixels;
                    event_x = round(events{i}{2}/100 * (screen_width - diag_max - border)) + diag_max/2 + border/2;
                    overlay_x = (1:resized_rotated_x)+double(round(round(x_padding/2)+event_x-resized_rotated_x/2));
                    event_y = round(events{i}{3}/100 * (screen_height - diag_max - border)) + diag_max/2 + border/2;
                    overlay_y = (1:resized_rotated_y)+double(round(screen_height+round(y_padding/2)-event_y-resized_rotated_y/2));
                    
                    alpha_mask(overlay_y,overlay_x) = resized_rotated_alpha;
                    image_overlay(overlay_y,overlay_x,:) = resized_rotated_event;
                    break; % found the image
                end
            end
        end
        y_range = (1:screen_height)+round(y_padding/2);
        x_range = (1:screen_width)+round(x_padding/2);
        image_overlay = image_overlay(y_range,x_range,:);
        alpha_mask = alpha_mask(y_range,x_range);
    end
end