%Code to track motion of particles
clearvars
clc

warning off

%file = '../data/002_AAC__20230921_active10_400nm_probtip_014_crop.nd2';

file = '../data/002_AAC__20230921_Passive_400nm_probtip_001.nd2';

reader = BioformatsImage(file);

trackerBG = LAPLinker;

%Initialize struct to hold particle information
particlePos = struct;

vid = VideoWriter('../processed/20231122_passive.avi');
vid.FrameRate = 5;
open(vid)

for iT = 1:reader.sizeT

    I = getPlane(reader, 1, 1, iT);

    %Register image
    if iT == 1

        refImg = I;

    else

        pxShift = xcorrreg(refImg, I);

        I = circshift(I, pxShift);
        
        refImg = I;

    end

    %Segment particle
    g1 = imgaussfilt(I, 3);
    g2 = imgaussfilt(I, 6);
    Idiff = g1 - g2;
    
    spotMask = Idiff > 2000;
    %imshowpair(Idiff, spotMask);

    %Get position, then fit a Guassian to better localize
    spotData = regionprops(spotMask, 'Centroid');

    if numel(spotData) >= 1

        %Crop the region around the spot
        xSpot = round(spotData(1).Centroid(1));
        ySpot = round(spotData(1).Centroid(2));
        spotImgCrop = I((ySpot - 15):(ySpot + 15), (xSpot - 15):(xSpot + 15));
        %imshow(spotImgCrop, [])

        %Coordinates for fitting
        xFit = 1:size(spotImgCrop,2);
        xFit = xFit - median(xFit);

        yFit = 1:size(spotImgCrop, 1);
        yFit = yFit - median(yFit);

        [xFit,  yFit] = meshgrid(xFit, yFit);

        %opts = optimset('display','off');

        fittedParams = lsqnonlin(@(params) fngauss(params, xFit, yFit) - double(spotImgCrop), ...
            [double(max(spotImgCrop(:))), 0, 0, 3, 3, mean(double(spotImgCrop(:)))]);

        % surf(xFit, yFit, spotImgCrop)
        % hold on
        % plot3(xFit, yFit, fngauss(fittedParams, xFit, yFit), 'o')
        % hold off

        spotX = xSpot + fittedParams(2);
        spotY = ySpot + fittedParams(3);

    else

        spotX = NaN;
        spotY = NaN;

    end


    %Collect data
    particlePos(iT).Centroid = [spotX, spotY];

    % %Check results
    % imshow(I, [])
    % hold on
    % plot(spotX, spotY, 'o')
    % hold off

    %Segment background beads
    Ibg = imgaussfilt(I, 2);
    Ibg = imopen(Ibg, strel('disk', 5));

    [centers, radii] = imfindcircles(Ibg, [5 20], 'Sensitivity', 0.88);

    %Remove the bright spot
    distToParticle = sqrt(sum((centers - [spotX, spotY]).^2, 2));

    idxParticle = find(distToParticle <= 2);

    centers(idxParticle, :) = [];
    radii(idxParticle, :) = [];

    %Compile the data into a struct for tracking
    for iBG = 1:size(centers, 1)
        bgData(iBG).Centroid = centers(iBG, :);
    end

    trackerBG = assignToTrack(trackerBG, iT, bgData);

    %Make the output image
    Iout = double(I);
    Iout = (Iout - min(Iout(:)))/(20000 - min(Iout(:)));
    Iout(Iout > 1) = 1;
    Iout(Iout < 0) = 0;
    %imshow(Iout)

    if ~isempty(centers)
        Iout = insertShape(Iout, 'circle', [centers(:, 1), centers(:, 2), ones(size(centers, 1), 1) * 2]);
    end

    if ~isnan(spotX)
        Iout = insertShape(Iout, 'filledcircle', [spotX, spotY, 8], 'color', 'magenta');
    end

    Iout = Iout ./ (max(Iout(:)));
    
    %imshow(Iout)
    writeVideo(vid, Iout)

end
close(vid)

save('../processed/20231122_passive.mat', 'trackerBG', 'particlePos')