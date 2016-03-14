
%This code performs the segmentations for all the cells in the fov for all
%the time steps
%Author: Tan H. Nguyen
%University of Illinois at Urbana-Champaign
clc;
clear all;
close all;
outdir = 'Z:\Hela_cell_time_laps_Feb_16th_2016\';
f=0:0;
t=88:88;
chh=0;
ii=0;
r=0:24;
z=0:0;
c =0:24;
%File name for the Halo-removed image
fout_slim_hr_ns=@(odir,f,t,i,ch,c,r,z) sprintf('%s\\f%d_t%d_i%d_ch%d_c%d_r%d_z%d_HR_NS.tif',odir,f,t,i,ch,c,r,z); % FOV, TIME, Channel, Frame Number, PAT

%File name for the segmentation image
fout_slim_seg=@(odir,f,t,i,ch,c,r,z) sprintf('%s\\f%d_t%d_i%d_ch%d_c%d_r%d_z%d_SEG.tif',odir,f,t,i,ch,c,r,z); % FOV, TIME, Channel, Frame Number, PAT

 for ff=f
      for rr=r   
           for cc=c 
               for tt = t
                    for zz = z
                           fnsname = fout_slim_hr_ns(outdir,ff,tt,ii,chh,cc,rr,zz);
                           fsegname = fout_slim_seg(outdir,ff,tt,ii,chh,cc,rr,zz);
                           im = imread(fnsname);
                           [bw_dil,ncells,bound_map]=cell_seg(im);
                           figure(1)
                           imagesc(im);
                           hold on;
                           for cellidx = 1:ncells
                               boundary = bound_map{cellidx};
                               plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
                           end
                           bw_dil = cast(bw_dil,'uint16');
                           %Write the binary mask for the cell
                           writeTIFF(bw_dil,fsegname,'uint16');
                           
                    end
               end
           end
      end
 end
 
