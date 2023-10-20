function [thresh_image] = ThresholdImageCalculator()
%function [thresh_image] = ThresholdImageCalculator()
%
%This function uses a set of foreground images to create an image to use as
%the threshold image. All foreground images must have the same size and
%align with one another. Once you get the image you will have to add back
%some pixels to make it agree with the camera width.

    %first load the foreground images
    [FileName, DataPath] = uigetfile('*.*','Select the images you would like to use to calculate the background', 'MultiSelect', 'on');

    %second cycle through the foreground images and store the maximum pixel
    %value
    thresh_image = imread([DataPath FileName{1}]);
    min_pix=10000;
    for j = 2:length(FileName)
        temp=imread([DataPath FileName{j}]);
        thresh_image = max(thresh_image, temp);
        min_pix=min(min_pix,min(min(temp)));
    end

    %next dilate the maximum pixel value over 10 pixels
    pixel_dilate=30;
    SE=strel('square',pixel_dilate); %structure element
    thresh_image_dilate = imdilate(thresh_image,SE);
      
    %then blur the image with a gaussian filter
    h = fspecial('gaussian',100,100); % params are filter type, size, sigma
    thresh_image_blur=imfilter(thresh_image_dilate,h);
    
    %finally normalize the image
    thresh_image_norm=double(thresh_image_blur-min_pix)/double(max(max(thresh_image_blur))-min_pix);

    %output the image
    thresh_image=thresh_image_norm;
end