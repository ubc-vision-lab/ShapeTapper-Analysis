function event_info = calculate_event_info(screen_w, screen_h, screen_dpi, safety, shape, img_dim, border_size_inches, event_x, event_y, event_rotation)
    new_screenH = screen_h * safety / 100;
    img_diag = norm(img_dim,2);
    toScale = img_diag / new_screenH; %toScale back up
    img_diag = img_diag / toScale;

    border = screen_dpi * border_size_inches; %0.375 inches = ~ 1cm (so 5mm border around to screen)
    border = 0; %screen_dpi/(screen_h/5/2); % Camera.getOrthogonalSize = 5 (from Unity)
    
    margin = screen_w - img_diag - border;
    toShift = event_x * margin / 100;
%     shiftAmountX = toShift - (margin / 2);
%     shiftAmountX = round(toShift - (margin / 2));
    shiftAmountX = toShift + (img_diag + border)/2; % gives us the absolute space
%     shiftAmountX = round(toShift + (img_diag + border)/2); % gives us the absolute space

    margin = screen_h - img_diag - border;
    toShift = event_y * margin / 100;
%     shiftAmountY = toShift - (margin / 2);
%     shiftAmountY = round(toShift - (margin / 2));
    shiftAmountY = toShift + (img_diag + border)/2;
%     shiftAmountY = round(toShift + (img_diag + border)/2);

    plusX = img_dim(1)/2;
    plusY = img_dim(2)/2;
    event_info = {toScale, -1 * shiftAmountX, -1 * shiftAmountY, 360-event_rotation, shape, plusX, plusY};
end