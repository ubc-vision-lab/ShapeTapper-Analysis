%%%%%%%%%%%%%%%%%%%%5
% Outputs the position and rotation for one image in a trial
% Inputs:
%       Given by Unity (demographics):
%       window_w, window_h, screen_dpi: Usually given by the demographics
%           file, in pixels
%       scaling: percentage of screen height to
%       shape: name of shape we want
%       img_dim: pixel resolution of image (x,y)
%       border_size_inches: extra margin in the trial, in inches
%       event_x, event_y: x,y position in percentage, relative to allowable space for
%       image (it's complicated)
%       event_rotation: degree rotation counterclockwise in Unity
%
% Output: Cell array containing the following, in order:
%       to_scale: The scaling factor
%       img_unity_x: The negative of the image shifting in x, in pixel space
%       img_unity_y: The negative of the image shifting in y, in pixel space
%       360 - event_rotation: CW rotation of the image (Unity is CCW)
%       shape: name of the shape
%       plusX, plusY: half the x and y resolutions of the image, so that
%       later transformations can use it for image-relative coordinates
%       (rather than from the centre)

function event_info = CalculateEventInfo(window_width, window_height, screen_dpi, scaling, shape, img_dim, border_size_inches, event_x, event_y, event_rotation)
    new_img_diag = window_height * scaling / 100; % what was the diagonal in Unity
    to_scale = norm(img_dim,2) / new_img_diag; % scaling factor so that original image diagonal will match new image diagonal

    border = screen_dpi * border_size_inches; %0.375 inches = ~ 1cm (so 5mm border around to screen)
    
    allowable_x = window_width - new_img_diag - border; % space where image is allowed, x coordinate
    img_unity_x = event_x * allowable_x / 100 + (new_img_diag + border)/2;

    allowable_y = window_height - new_img_diag - border;
    img_unity_y = event_y * allowable_y / 100 + (new_img_diag + border)/2;

    % image coordinate correction
    plusX = img_dim(1)/2;
    plusY = img_dim(2)/2;
    
    % Unity rotates CCW, MATLAB rotates, CW. Subtract rotation from 360 to
    % get the matching direction
    event_info = {to_scale, img_unity_x, img_unity_y, 360-event_rotation, shape, plusX, plusY};
end