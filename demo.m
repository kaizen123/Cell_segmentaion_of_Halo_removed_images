%This script demonstrate the cell segmentation results
%Author: Tan H. Nguyen
%University of Illinois at Urbana-Champaign
function demo
    clc;
    clear all;
    close all;
    im = cast(imread('D:\Hela_cell_time_laps_Feb_16th_2016\f0_t0_i0_ch0_c6_r13_z0_HR_NS.tif'),'single');

    tic;
    [~,ncells,bound_map]=cell_seg(im);
    toc;
    figure(2);
    imagesc(im);colorbar;
    hold on;
    for cellidx = 1:ncells
       boundary = bound_map{cellidx};
       plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
    end

end
