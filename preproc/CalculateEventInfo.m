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

function [A] = CalculateEventInfo(scrnWidth, scrnHeight, dpi, unityScreenMargin, ...
    rawImgsWidthsAndHeights, scaling, rotation, XPosRawImgsCntrPrcnt,...
    YPosRawImgsCntrPrcnt)

%compute the diagonal length of the image presented in unity (in pixels). Note that this
%is a proportion of the screen height   
unityImgDiags = transpose(scrnHeight*(scaling/100));

unityMarginPixels = dpi*unityScreenMargin; %0.375 inches = ~ 1cm (so 5mm border around to screen)

%compute the scaling factor to convert the raw image to unity size
scaleUnty2RawImgs = sqrt(sum(rawImgsWidthsAndHeights.^2,2))./unityImgDiags; 

%compute the total margin along the x- and y-coordinates (screen width and height) in unity pixels
totXMargin = scrnWidth - unityImgDiags - unityMarginPixels;
totYMargin = scrnHeight - unityImgDiags - unityMarginPixels;

%compute the x and y positions of the centre of the images in unity space (in pixels)
XPosCntrUntyImgs = ((totXMargin.*transpose(XPosRawImgsCntrPrcnt/100)) + ...
    ((unityImgDiags + unityMarginPixels)/2));
YPosCntrUntyImgs = ((totYMargin.*transpose(YPosRawImgsCntrPrcnt/100)) + ...
    ((unityImgDiags + unityMarginPixels)/2));

%image coordinate correction to convert the image-centred TP coordinates back to a canvas
%based one.
plusX = rawImgsWidthsAndHeights(:,1)/2;
plusY = rawImgsWidthsAndHeights(:,2)/2;
    
% Unity rotates CCW, MATLAB rotates, CW. Subtract rotation from 360 to
% get the matching direction
A = [scaleUnty2RawImgs XPosCntrUntyImgs YPosCntrUntyImgs transpose(360-rotation) plusX plusY];
return;