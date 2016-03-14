
function [bw_dil,ncells,bound_map]=cell_seg(im)
    %Perform segmentation given the input image
    %Inputs:
    %   im: an input image (grayscale)
    %Outputs:
    %   bw_dil: a binary map of the input image
    %   k: number of cells currently in the field of view
    %   boundary: a structure containing the boundary of the segmented
    %   glands
     %Smooth out the image before moving on to avoid the noise
    display = 0;
    ho = fspecial('gaussian',[6 6],1);
    im = imfilter(im,ho,'same');
    %Generate the binary mask for the cells
    [~,threshold]=edge(im,'sobel');
    fudgeFactor = 0.4;%Smaller, pick more smaller change in the edge
    bw = edge(im,'sobel',threshold*fudgeFactor);
    %Dilate the image
    se90 = strel('line',4,90);
    se = strel('line',4,0);
    bw_dil = imdilate(bw,[se90 se]); %Dilate the image horizontally and vertically
    bw_dil = imfill(bw_dil,'holes');%Fill out the holes inside the cells
    minarea = 2000;

    %Make sure the background is nicer
    bg = 1-bw_dil;
  
    bg_dil = imclose(bg,strel('disk',5));
    bg_dil = bwareaopen(bg_dil,50);
    bw_dil = bw_dil.*(1-bg_dil);
    bw_dil = bwareaopen(bw_dil,minarea);%Make sure we get rid of too small cells
    L = bwlabel(bw_dil,4);
   
    if (display)
        subplot(222);
        imagesc(bw_dil);colormap gray;

        subplot(223);
        imagesc(bg_dil);title('Background mask');
        
        subplot(224);
        imagesc(L);colormap jet;
   end
    bound_map = bwboundaries(bw_dil,4);
    ncells =size(bound_map,1);
    
end