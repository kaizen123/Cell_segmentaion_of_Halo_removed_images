%This script demonstrate the cell segmentation results
%Author: Tan H. Nguyen
%University of Illinois at Urbana-Champaign
function demo
    clc;
    clear all;
    close all;
    im = cast(imread('E:\Data_for_cell_segmentation\f0_t71_i0_ch0_c10_r16_z0_HR_NS.tif'),'single');

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
