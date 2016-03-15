%This script demonstrate the cell segmentation results
%Author: Tan H. Nguyen
%University of Illinois at Urbana-Champaign
function demo
    clc;
    clear all;
    close all;
    im = cast(imread('C:\Users\QLI\Desktop\Full_fov_segmentation_images\f0_t0_i0_ch0_c12_r14_z0_HR_NS.tif'),'single');

    tic;
    [~,ncells,bound_map]=cell_seg(im);
    toc;
    figure(2);
    imagesc(im);
    hold on;
    for cellidx = 1:ncells
       boundary = bound_map{cellidx};
       plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
    end

end
