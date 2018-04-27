function [A] = ProcessSTFiles(readOrWrite,storeImages)

switch nargin
    case 0
        readOrWrite = 'r';
        storeImages = 'n';
    case 1
        storeImages = 'n';
end
        
%FIXED VARIABLES

%image file extension
IMGFEXT = '.png';

%unity screen margin (in inches)
UNTYSCRNMARGIN = 0;

%length of subject ID name in characters
SIDNMLNGTH=4;
BLCKNMFNINDX1 = 6;
BLCKNMFNINDX2 = 7;

%Configuration file format specifiers handle three cases (1 image per trial,
%2 images per trial, or 3 images per trial). Function requires the first
%line of the configuration file to be analyzed to determine this. The final
%four fields are all 0s (one image per trial), a pair of zeros (two images
%per trial), or all non-0s (three images per trial). It is assumed when two
%images are presented per trial, the first 2 elements of the final 4 will
%be non-zero.
FORMATCNFG = ['%f%f%*f%*f%*f%*f%*f' repmat('%*s',[1 6]) '%*f%f%f'...
    repmat('%f%s%f%f%*f%*f%*f%*f',[1 3]) '%*f%f%f%f%f'];

%The fields, columns, in the results files are as follows:
%1: trial # within the block; 2: Badflag (0=No Error; 2=Error); 3: RT; 4: x touchpoint in
%pixels; 5: y touchpoint in pixels.

%format specifier for the reuslts files
FORMATRSLTS = '%f%f%f%f%f';

    %merge the configuration data and results file data into one cell data structure to
    %store to disk. Structure uses rows as fields for 'fprintf' as follows
    %1: Trial #; 2: Block #; 3:Trial # within the block; 4:BadFlag; 5: RT;
    %6: x-touchpoint (pixels); 7: y-touchpoint (pixels); 8: Event 1 Image Name; 9: Onset
    %Delay; 10: E1 scaling factor (%); 11: E1 Rotation; 12: E1 x-centre of image (% of screen);
    %13: E1 y-centre of image (% of screen); 14: Event 2 Image Name; 15: E2 Onset Delay; 
    %16: E2 scaling factor (%); 17: E2 Rotation; 18: E2 x-centre of image (% of screen);
    %19: E2 y-centre of image (% of screen); 20: E3 Image Name; 21: E3 Onset
    %Delay; 22: E3 scaling factor (%); 23: E3 Rotation; 24: E3 x-centre of image (% of screen);
    %25: E3 y-centre of image (% of screen); 
    
%Header for the stored data
HDR = ['Trial #\tBlock #\tTrial # w/n block\tBad Flag\tRT\tX-TP (pixels, 0,0 is bottom left)'...
        '\tY-TP (pixels, 0,0 is bottom left)\tEvent 1 Image Name\tOnset Delay\tscaling factor (%%)'...
        '\tE1 Rotation (degrees CCW)\tE1 x-centre of image (%% of screen width)\t'...
        'E1 y-centre of image (%% of screen height)\tE2 Image Name\tE2 Onset Delay'...
        '\tE2 scaling factor (%%)\tE2 Rotation\tE2 x-centre of image (%%)\t'...
        'E2 y-centre of image (%%)\tE3 Image Name\tE3 Onset Delay\tE3 scaling factor (%%)\t'...
        'E3 Rotation\tE3 x-centre of image (%%)\tE3 y-centre of image (%%)\t'...
        'Transformed X-TP (pixels, 0,0 is bottom left)\tTransformed Y-TP (pixels, 0,0 '...
        'is bottom left)\tClosest Image\tOutside Canvas?\r\n'];

%format specifier for the merged results and confirugation files that are
%stored to disk
FORMATSAVEDDATA = ['%f' repmat('\t%f',[1 6]) repmat('\t%s\t%f\t%f\t%f\t%f\t%f',[1,3]) ...
    '\t%f\t%f\t%s\t%f\r\n'];

%use whichever directory you want to work from
rootDirPath = uigetdir();
cd(rootDirPath);

%select the demographic files
[demoFNmLst demFDir] = uigetfile('*demographic.txt','MultiSelect','On');

%determine the number of demographic files and ensure file names are stored in a cell array
if iscell(demoFNmLst )
    numDemFls = length(demoFNmLst);
else
    numDemFls = 1;
    demoFNmLst = {demoFNmLst};
end

demogFDirPath = fullfile(rootDirPath,'Demographics');
rawDataFDirPath = fullfile(rootDirPath,'Raw_Results_Files');
cnfgFDirPath = fullfile(rootDirPath,'Configs');
rawImgFDirPath = fullfile(rootDirPath,'Images');
trialImgFDirPath = 'Trial_Image_Results';

%navigate to the folder that holds the demographic files.
cd(demogFDirPath)

%extract and store the image files
cd(rawImgFDirPath)
rawImgFLst = ls(['*' IMGFEXT]);

%loop through each image file name in the list to load and store the image
%name and its height and width dimensions into a dictionary
rawImgNames = cell(1,size(rawImgFLst,1));
rawImgWidthsAndHeights = rawImgNames;
rawImgs = rawImgNames;
rawImgsAlpha = rawImgNames;

for i = 1:size(rawImgFLst,1)
    
    currImgName = strtrim(rawImgFLst(i,:));
    
    %store the image info
    [rawImgs{i},~,rawImgsAlpha{i}] = imread(currImgName,'BackgroundColor','none');
    
    %store the height and width of the image (in pixels)
    [y,x,~] = size(imread(currImgName));
    rawImgWidthsAndHeights{i} = [x y];
    
    %store the image name
    rawImgNames{i}=currImgName(1:end-4);
end

%create the dictionary
rawImgSizeDict = containers.Map(rawImgNames,rawImgWidthsAndHeights);
rawImgsRGBDict = containers.Map(rawImgNames,rawImgs);
rawImgsAlphaDict = containers.Map(rawImgNames,rawImgsAlpha);

%loop through each demographic file to extract the experimental parameters:
%1: user ID; 2: the configuration file; 3: gender; 4: handedness; 5: age;
%6: screen width resolution in pixels; 7: screen height resolution in pixels;
%8: DPI; 9: configuration file name

for i = 1:numDemFls
    
    %navigate to the demographic file directory
    cd(demogFDirPath)
    
    %open the current demographic file
    fid = fopen(demoFNmLst{i}, 'r');
    demoData = textscan(fid,'%s','Delimiter',',','HeaderLines',1);
    fclose(fid);
    
    %extract and store the subject ID
    sID = demoFNmLst{i}(1:SIDNMLNGTH);
    
    %extract and store the experimental paramaters
    scrnWdthPxls = str2double(demoData{1}{5,:}); %screen width of device in pixels
    scrnHghtPxls = str2double(demoData{1}{6,:}); %screen height of device in pixels
    dpi = str2double(demoData{1}{7,:}); %screen DPI of the device
    cnfgFlNm = demoData{1}{8,:}; %config file name which includes the file extension
    
    %navigate to the raw data folder
    cd(rawDataFDirPath);
    
    %extract and store a list of the names of the results files
    ssRsltsFNmLst = ls(['*' sID '*.txt']);
    
    %loop through each result file with the current demographic user ID name and extract
    %the content of the file (the results data) and the block numbers (each block was
    %stored as a separate file) from the file name
    resultsData=[];
    for j=1:size(ssRsltsFNmLst,1)
        
        fid = fopen(ssRsltsFNmLst(j,:),'r');
        currResults = cell2mat(textscan(fid,FORMATRSLTS,'Delimiter',',','HeaderLines',1));
        fclose(fid);
        
        %if the block number is less than 10, it will have a '.' in the second element
        if strcmp(ssRsltsFNmLst(j,BLCKNMFNINDX2),'.')
            blckNum = str2num(ssRsltsFNmLst(j,BLCKNMFNINDX1));
        else
            blckNum = str2num(ssRsltsFNmLst(j,BLCKNMFNINDX1:BLCKNMFNINDX2));
        end
        
        %prepend the block number to 'currResults'
        currResults = [blckNum*ones(size(currResults,1),1) currResults];
        
        %update 'rsltsData'
        resultsData=[resultsData; currResults];
    end
    
    %store the trial total for the set of results files associated with the current
    %demographic file
    trialTot = size(resultsData,1);
    
    %extract and store the configuration file information
    cd(cnfgFDirPath);
    fid = fopen(cnfgFlNm,'r');
    cnfgData = textscan(fid,FORMATCNFG,'Delimiter',',');
    fclose(fid);
    
    %merge the configuration data and results file data into one cell data structure to
    %store to disk. Structure uses rows as fields for 'fprintf' as follows
    %1: Trial #; 2: Block #; 3:Trial # within the block; 4:BadFlag; 5: RT;
    %6: x-touchpoint (pixels); 7: y-touchpoint (pixels); 8: Event 1 Image Name; 9: Onset
    %Delay; 10: E1 scaling factor (%); 11: E1 Rotation; 12: E1 x-centre of image (% of screen);
    %13: E1 y-centre of image (% of screen); 14: Event 2 Image Name; 15: E2 Onset Delay; 
    %16: E2 scaling factor (%); 17: E2 Rotation; 18: E2 x-centre of image (% of screen);
    %19: E2 y-centre of image (% of screen); 20: E3 Image Name; 21: E3 Onset
    %Delay; 22: E3 scaling factor (%); 23: E3 Rotation; 24: E3 x-centre of image (% of screen);
    %25: E3 y-centre of image (% of screen); 26: transformed X-TP (pixels); 27: transformed
    %Y-TP (pixels); 28: Closest Image (1-3); 29: Within Canvas.
    
    %The columns (fields) in 'resultsData' (the the results files content) are as follows:
    %1: block #; 2: trial # within the block; 3: Badflag (0=No Error; 2=Error); 4: RT;
    %5: x touchpoint (pixels); 6: y touchpoint (pixels).
 
    %The fields, columns, in 'cnfgData' (the content of the configuration files) are as follows:
    %1: Block #; 2: Trial # within the block; 3: Event 1 (E1) x-centre of image (% of screen);
    %4: E1 y-centre of image (% of screen); 5: E1 Onset Delay; 6: EI Image Name;
    %7: E1 Rotation; 8: E1 scaling factor (%); 9: E2 Onset Delay; 10: E2 Image Name;
    %11: E2 Rotation; 12: E2 scaling factor (%); 13: E3 Onset Delay; 14: E3 Image Name;
    %15: E3 Rotation; 16: E3 scaling factor; 17: E2 x-centre of image (% of screen);
    %18: E2 y-centre of image (% of screen); 19: E3 x-centre of the image (% of screen);
    %20: E3 y-centre of image (% of screen);
    
    %merge the configuration file and the block
    mrgdData = cell(29,trialTot);
    
    %a matrix comprised of the heights and widths of the raw images for
    %each event. A nx2x3 3D matrix with rows=trials, columns (width, height), and
    %pages reflect events 1,2 and 3.
    rawImgsWidthsAndHeights = nan(size(resultsData,1),2,3);
    
    cnfgBlckPerTrial=cnfgData{1}(:); 
    currRsltsBlck = 0;
    
    %loop through the total number of trials (cols of 'mrgdData')
    for j = 1:trialTot
        
        %store the next block number from the 'resultsData' 
        nxtRsltsBlck = resultsData(j,1);
        
        %if we've moved on to a different block number in the results data
        %set, start a new row index for 'cnfgData'
        if currRsltsBlck~=nxtRsltsBlck
            
            %the row index corresponding to the first entry of the next
            %block number
            iniCnfgRow = find(cnfgBlckPerTrial==nxtRsltsBlck,1);
 
            %update the current block number
            currRsltsBlck=nxtRsltsBlck;
        end
        
        %update the 'cnfgData' index ('cnfgRow')
        cnfgRow = iniCnfgRow + resultsData(j,2) - 1;
   
        %store all of the variables to 'mrgdData'
        mrgdData{1,j} = j; %overall trial counter
        mrgdData{2,j} = resultsData(j,1); %Block #
        mrgdData{3,j} = resultsData(j,2); %Trial # within the block
        mrgdData{4,j} = resultsData(j,3); %badflag
        mrgdData{5,j} = resultsData(j,4); %RT
        mrgdData{6,j} = resultsData(j,5); %x-touchpoint (pixels)
        mrgdData{7,j} = resultsData(j,6); %y-touchpoint (pixels)
        mrgdData{8,j} = cnfgData{6}{cnfgRow}; %E1 Image Name
        mrgdData{9,j} = cnfgData{5}(cnfgRow); %E1 Onset Delay
        mrgdData{10,j} = cnfgData{8}(cnfgRow); %E1 scaling factor (%)
        mrgdData{11,j} = cnfgData{7}(cnfgRow); %E1 rotation
        mrgdData{12,j} = cnfgData{3}(cnfgRow); %E1 x-centre of image (% of screen)
        mrgdData{13,j} = cnfgData{4}(cnfgRow); %E1 y-centre of image (% of screen)
        rawImgsWidthsAndHeights(j,:,1) = rawImgSizeDict(cnfgData{6}{cnfgRow});

        %if there is no event 2 (check image name, 'NaN' means no Event)
        if isempty(cnfgData{10}{cnfgRow})
            mrgdData{14,j} = NaN;
            mrgdData{15,j} = NaN;
            mrgdData{16,j} = NaN;
            mrgdData{17,j} = NaN;
            mrgdData{18,j} = NaN;
            mrgdData{19,j} = NaN;
            rawImgsWidthsAndHeights(j,:,2) = [NaN NaN];
        else
            mrgdData{14,j} = cnfgData{10}{cnfgRow}; %E2 Image Name
            mrgdData{15,j} = cnfgData{9}(cnfgRow); %E2 onset delay
            mrgdData{16,j} = cnfgData{12}(cnfgRow); %E2 scaling factor (%)
            mrgdData{17,j} = cnfgData{11}(cnfgRow); %E2 rotation
            mrgdData{18,j} = cnfgData{17}(cnfgRow); %E2 x-centre of image (% of screen)
            mrgdData{19,j} = cnfgData{18}(cnfgRow); %E2 y-centre of image (% of screen)
            rawImgsWidthsAndHeights(j,:,2) = rawImgSizeDict(cnfgData{10}{cnfgRow});
        end
        
        %if there is no event 3 (check image name, 'NaN' means no Event)
        if isempty(cnfgData{14}{cnfgRow})
            mrgdData{20,j} = NaN;
            mrgdData{21,j} = NaN;
            mrgdData{22,j} = NaN;
            mrgdData{23,j} = NaN;
            mrgdData{24,j} = NaN;
            mrgdData{25,j} = NaN;
            rawImgsWidthsAndHeights(j,:,3) = [NaN NaN];
        else
            mrgdData{20,j} = cnfgData{14}{cnfgRow}; %E3 Image Name
            mrgdData{21,j} = cnfgData{13}(cnfgRow); %E3 onset delay
            mrgdData{22,j} = cnfgData{16}(cnfgRow); %E3 scaling factor (%)
            mrgdData{23,j} = cnfgData{15}(cnfgRow); %E3 rotation
            mrgdData{24,j} = cnfgData{19}(cnfgRow); %E3 x-centre of image (% of screen)
            mrgdData{25,j} = cnfgData{20}(cnfgRow); %E3 y-centre of image (% of screen)
            rawImgsWidthsAndHeights(j,:,3) = rawImgSizeDict(cnfgData{14}{cnfgRow});
        end       
    end

    %compute initial tranformation parameters for each event. 'eventInfo' rows reflect trials
    %and coliumns reflect CalculateEventInfo output variables 'scaleRawImgs2Unty',
    %XPosCntrUntyImgs, YPosCntrUntyImgs, Matlab rotation (degrees), plusX, plusY
    eventInfo = nan(trialTot,6,3);
    eventInfo(:,:,1) = CalculateEventInfo(scrnWdthPxls,scrnHghtPxls,dpi,UNTYSCRNMARGIN,...
        rawImgsWidthsAndHeights(:,:,1),[mrgdData{10,:}],[mrgdData{11,:}],...
        [mrgdData{12,:}],[mrgdData{13,:}]);
    eventInfo(:,:,2) = CalculateEventInfo(scrnWdthPxls,scrnHghtPxls,dpi,UNTYSCRNMARGIN,...
        rawImgsWidthsAndHeights(:,:,2),[mrgdData{16,:}],[mrgdData{17,:}],...
        [mrgdData{18,:}],[mrgdData{19,:}]);
    eventInfo(:,:,3) = CalculateEventInfo(scrnWdthPxls,scrnHghtPxls,dpi,UNTYSCRNMARGIN,...
        rawImgsWidthsAndHeights(:,:,3),[mrgdData{22,:}],[mrgdData{23,:}],...
        [mrgdData{24,:}],[mrgdData{25,:}]);
    
    %the distance between the touchpoint and the centre(s) of the image(s)
    scaledDistTPs2CentreUntyImg = nan(trialTot,3);
    untyImgRadii = nan(trialTot,3);
    
    %initialize 'scaledTPs' which will store the TPs scaled to unity space
    scaledTPs = nan(trialTot,2,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   TRANSFORM THE DATA TO AN CENTRE-OF-THE-IMAGE IN UNITY AND MATLAB SPACE   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %store the touchpoints and centre positions for the unity images
    rawTouchPoints = [transpose([mrgdData{6,:}]) transpose([mrgdData{7,:}])]; %x (width) then y (height)
    
    %loop through each unity image events 1 to 3
    for j=1:3
        
        %centre positions of the unity images, x (width) then y (height)
        untyImgsCentrePos = [eventInfo(:,2,j) eventInfo(:,3,j)]; 
    
        %scale the unity-image-centred touchpoints to the raw image space
        scaledTPs(:,:,j) = (rawTouchPoints - untyImgsCentrePos).*[eventInfo(:,1,j) eventInfo(:,1,j)];   
    
        %compute the distance between the touchpoints and the centre of the unity
        %images
        scaledDistTPs2CentreUntyImg(:,j) = sqrt(sum(scaledTPs(:,:,j).^2,2));
        
        %the radii of the unity images
        untyImgRadii(:,j) = sqrt(sum([eventInfo(:,5,j) eventInfo(:,6,j)].^2,2));
    end
    
    %determine which image ('event') the touchpoint is closest to
    [~, closestImg] = min(scaledDistTPs2CentreUntyImg,[],2);
 
    %extract the angles, plus Xs and plus Ys associated with the closest TPs
    lnInx = sub2ind([trialTot 3],transpose(1:trialTot),closestImg);
    allAngles = squeeze(eventInfo(:,4,:));
    clstImgAngles = allAngles(lnInx); 
    
    %compute the rotation angles for the scaled TPs and store them in a rotation matrix.
    %Each page is a trial.
    R2 = nan(2,2,trialTot);
    R2(1,1,:) = cosd(clstImgAngles);
    R2(1,2,:) = sind(clstImgAngles);
    R2(2,1,:) = -sind(clstImgAngles);
    R2(2,2,:) = cosd(clstImgAngles);
    
    %initialize 'trnsfrmdTPs' extract and store the scaled TPs that are closest to the image
    trnsfrmdTPs = nan(trialTot,2);
    
    %update 'mrgdData' with the tansformed touchpoints, image assignments, and an
    %evaulation of whether or not the touchpoint fell within the image canvas
    for j=1:trialTot
        trnsfrmdTPs(j,:) = scaledTPs(j,:,closestImg(j))*R2(:,:,j);
        
        %add half the width of the target image to the X-TP
        trnsfrmdTPs(j,1) = trnsfrmdTPs(j,1) + eventInfo(j,5,closestImg(j));
        
        %add half the heighth of the target image  to the Y-TP
        trnsfrmdTPs(j,2) = trnsfrmdTPs(j,2) + eventInfo(j,6,closestImg(j));
        
        mrgdData{26,j} = trnsfrmdTPs(j,1); %transformed X-touchpoint
        mrgdData{27,j} = trnsfrmdTPs(j,2); %transformed Y-touchpoint
        
        %if the TP has a closest image
        if ~isnan(closestImg(j))
            
            switch closestImg(j)
                case 1 %event (image) 1
                    mrgdData{28,j} = mrgdData{8,j}; %closest image to the touchpoint is E1
                case 2 %event (image) 2
                    mrgdData{28,j} = mrgdData{14,j}; %closest image to the touchpoint is E2
                case 3 %event (image) 2
                    mrgdData{28,j} = mrgdData{20,j}; %closest image to the touchpoint is E3
            end
            
            %if the scaled touchpoint distance is > the radius of the unity
            %image canvas, it is outside of the canvas
            if scaledDistTPs2CentreUntyImg(j,closestImg(j)) > untyImgRadii(closestImg(j))
                mrgdData{29,j} = 1;
            else
                mrgdData{29,j} = 0;
            end

        else
            mrgdData{28,j} = NaN;
            mrgdData{29,j} = NaN;
        end
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   SAVE ANALYSIS TO DISK (IF NEED BE)  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
        %if storing the merged data set is desirable
    if strcmp(readOrWrite,'w')
        
       %store the merged data file to disk
        cd(rootDirPath);
        
        %store the merged data set
        fid = fopen(['Results_' cnfgFlNm(1:end-4) '_' sID '.txt'],'w');

        %store the header, and then the remainder of the data
        fprintf(fid,HDR);
        fprintf(fid,FORMATSAVEDDATA,mrgdData{:,:});
        fclose(fid);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   STORE IMAGES PER TRIAL (IF NEED BE) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(storeImages,'y')
        
        %return to the root directory
        cd(rootDirPath);
        
        %check if the images per trial folder exists. If not, create one.
        if exist(trialImgFDirPath,'dir')==7
        else
            mkdir(trialImgFDirPath);     
        end

        %navigate to the images per trial folder
        cd(trialImgFDirPath); 

        %create a folder name specific to the current demographic file and its associated data
        %sets
        fldrNm = [cnfgFlNm(1:end-4) '_' sID];

        %check if a folder specific to the demographics file currently being procesesd exists.
        %If not, create one.
        if exist(fldrNm,'dir')==7
        else
            mkdir(fldrNm);     
        end

        %navigate to the folder
        cd(fldrNm);
        
        %the canonical images when rotated to match their presentation in unity can
        %sometimes end up partially out-of-bounds. Add padding to the
        %canvas.
        scrnPaddingPxls = 300;
        
        %create a touchpoint image and colour it (magenta)
        lengthGraphicTPPxls = 10;
        graphicTPRGB = 255*ones(lengthGraphicTPPxls,lengthGraphicTPPxls,3,'uint8');
        graphicTPRGB(:,:,2) = 0;
        graphicTPAlpha = 255*ones(lengthGraphicTPPxls,lengthGraphicTPPxls,'uint8');

        %re-create the unity window for each trial
        for j=1:trialTot
            
            %create a grey screen canvas with some additional padding beyond the unity window
            screenPlusPaddingRGB = ones(scrnHghtPxls+scrnPaddingPxls,scrnWdthPxls+scrnPaddingPxls, 3, 'uint8');
            screenPlusPaddingAlpha = zeros(scrnHghtPxls+scrnPaddingPxls,scrnWdthPxls+scrnPaddingPxls, 'uint8');
            
            %loop through each of three possible events (images) presented
            %in a given trial
            for k=1:3               
            	
                %if an image was presented in the current event, process it
                if ~isnan(mrgdData{8+((6*k)-6),j});
                
                    %retrieve the current raw image and its alpha
                	currRawImg = rawImgsRGBDict(mrgdData{8+((6*k)-6),j});
                	currRawImgAlpha = rawImgsAlphaDict(mrgdData{8+((6*k)-6),j});
                    
                    %scale the raw and alpha images to unity space using the inverse of 
                    %col 1 in 'eventInfo' (the scale factor from unity to raw) and rotate 
                    %the raw image using col 4 in 'eventInfo'
                    rawImgSzd2Unty = imresize(currRawImg,(1/eventInfo(j,1,k)));
                    rawImgAlphaSzd2Unty = imresize(currRawImgAlpha,(1/eventInfo(j,1,k)));
                    untyImg = imrotate(rawImgSzd2Unty,eventInfo(j,4,k));
                    untyImgAlpha = imrotate(rawImgAlphaSzd2Unty,eventInfo(j,4,k));
                    
                    %store the width and height of the unity image
                    [untyImgHghtPxls,untyImgWdthPxls,~] = size(untyImg); 
                    
                    %compute the location of the centre of the unity image as it would appear
                    %in matlab's coordinate system (where 0,0 is top left). This means
                    %that only the Y-centre needs to be converted to Matlab's system.
                    untyImgOnMtlbCvsCntrXPos = round(eventInfo(j,2,k)) + round(scrnPaddingPxls/2);
                    untyImgOnMtlbCvsCntrYPos = scrnHghtPxls - round(eventInfo(j,3,k)) + round(scrnPaddingPxls/2);
                    
                    untyImgOnMtlbCnvsWdthIndices = (1:untyImgWdthPxls)+untyImgOnMtlbCvsCntrXPos-round(untyImgWdthPxls/2);
                    untyImgOnMtlbCnvsHghtIndices = (1:untyImgHghtPxls)+untyImgOnMtlbCvsCntrYPos-round(untyImgHghtPxls/2);
                    
                    %add the unity image and the unity alpha image into the padded screen 
                    %space
                    screenPlusPaddingRGB(untyImgOnMtlbCnvsHghtIndices,untyImgOnMtlbCnvsWdthIndices,:) = untyImg;
                    screenPlusPaddingAlpha(untyImgOnMtlbCnvsHghtIndices,untyImgOnMtlbCnvsWdthIndices) = untyImgAlpha;               
                end
            end
            
            %if the touchpoint value is 0,0, do nothing because it's not valid
            if mrgdData{6,j}~=0 || mrgdData{7,j} ~= 0
                
                %compute the coordinates the touchpoint will occupy. Note that the unity
                %results files reflect the unity image coordinate system with (0,0) bottom
                %left. Matlab uses (0,0) as the top left. This requires changing the
                %unity y-touchpoint only
                
                matlabTPX = mrgdData{6,j} + round(scrnPaddingPxls/2);
                matlabTPY = (scrnHghtPxls + round(scrnPaddingPxls/2)) - mrgdData{7,j};
                matlabTPYIndices = (1:lengthGraphicTPPxls) + matlabTPY - round(lengthGraphicTPPxls/2);
                matlabTPXIndices = (1:lengthGraphicTPPxls) + matlabTPX - round(lengthGraphicTPPxls/2);
                
                %if j==117
                %    disp('stop!');
                %end
                %add the touchpoint to the canvas
                screenPlusPaddingRGB(matlabTPYIndices,matlabTPXIndices,:) = graphicTPRGB;
                screenPlusPaddingAlpha(matlabTPYIndices,matlabTPXIndices) = graphicTPAlpha;
            end
            
            %remove the padding
            untyScrnHghtIndcs = (1:scrnHghtPxls)+round(scrnPaddingPxls/2);
            untyScrnWdthIndcs = (1:scrnWdthPxls)+round(scrnPaddingPxls/2);
            screenRGB = screenPlusPaddingRGB(untyScrnHghtIndcs,untyScrnWdthIndcs,:);
            screenAlpha = screenPlusPaddingAlpha(untyScrnHghtIndcs,untyScrnWdthIndcs);
        
        %image filename
        imgFNm = ['Image_' sID '_' cnfgFlNm(1:end-4) '_'];
        
        %append the trial number to the end of the image filename
        if j < 10
            imgFNm = [imgFNm '00' num2str(j) IMGFEXT];
        elseif j > 9 && j < 100
            imgFNm = [imgFNm '0' num2str(j) IMGFEXT];
        else
            imgFNm = [imgFNm num2str(j) IMGFEXT];
        end
          
        %Save the image to disk
        imwrite(screenRGB,imgFNm,'png','Alpha',screenAlpha);  
        end
    end
    
    %store the merged results into the output variable 'A'.
    A{i} = mrgdData;
end
cd(rootDirPath);
return;

