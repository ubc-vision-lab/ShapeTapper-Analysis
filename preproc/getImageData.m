% Given a row of shapetapper configuration (a trial), generates the
% the image intended, at the input resolution
% input: trial configuration, resolution (x and y)
function [images, block, trial] = getImageData(trialConfig)
    %grab relevant columns
    % indices of columns. TODO: use a dictionary to index?
    block = str2double(trialConfig{1});
    trial = str2double(trialConfig{2});
    %not all images exist at any given point
    % in order: image name, x_position, y_position, rotation, safety
    image_indices = [[18,15,16,19,20];[26,42,43,27,28];[34,44,45,35,36]];
    images = cell(3,1);

    for i = 1:size(image_indices,1)
        if ~isempty(trialConfig{image_indices(i,1)})
            image_values = cell(1,5);
            image_values{1,1} = trialConfig{image_indices(i,1)};
            for j = 2:size(image_indices,2)
                image_values{1,j} = str2double(trialConfig{image_indices(i,j)});
            end
            images{i} = image_values;
        end
    end
    images = images(~cellfun('isempty',images));
end