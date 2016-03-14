
%This code performs the segmentations for all the cells in the fov for all
%the time steps
%Author: Tan H. Nguyen
%University of Illinois at Urbana-Champaign
clc;
clear all;
close all;
outdir = 'Z:\Hela_cell_time_laps_Feb_16th_2016\';
matlab_dir = strcat(outdir,'\mat_files\');
if (~exist(matlab_dir))
    mkdir(matlab_dir);
end
f=0:0;
t=[0:10 80:88];
chh=0;
ii=0;
r=0:24;
z=0:0;
c =0:24;
%File name for the Halo-removed image
fout_slim_hr_ns=@(odir,f,t,i,ch,c,r,z) sprintf('%s\\f%d_t%d_i%d_ch%d_c%d_r%d_z%d_HR_NS.tif',odir,f,t,i,ch,c,r,z); % FOV, TIME, Channel, Frame Number, PAT

%File name for the segmentation image
fout_slim_seg=@(odir,f,t,i,ch,c,r,z) sprintf('%s\\f%d_t%d_i%d_ch%d_c%d_r%d_z%d_SEG.tif',odir,f,t,i,ch,c,r,z); % FOV, TIME, Channel, Frame Number, PAT
fout_slim_overlaid=@(odir,f,t,i,ch,c,r,z) sprintf('%s\\f%d_t%d_i%d_ch%d_c%d_r%d_z%d_OVL.tif',odir,f,t,i,ch,c,r,z); % FOV, TIME, Channel, Frame Number, PAT
fout_bound_mat=@(odir,f,t,i,ch,c,r,z) sprintf('%s\\f%d_t%d_i%d_ch%d_c%d_r%d_z%d_bound.mat',odir,f,t,i,ch,c,r,z); % FOV, TIME, Channel, Frame Number, PAT

saveoverlaid = 0;
min_phaseval = -0.3;
max_phaseval = 2.5;
 for ff=f
      for rr=r   
           for cc=c 
               for tt = t
                    for zz = z
                           fnsname = fout_slim_hr_ns(outdir,ff,tt,ii,chh,cc,rr,zz);
                           fsegname = fout_slim_seg(outdir,ff,tt,ii,chh,cc,rr,zz);
                           fmatname = fout_bound_mat(matlab_dir,ff,tt,ii,chh,cc,rr,zz); %Name of the mat file to be saved
                           im = imread(fnsname);
                           [bw_dil,ncells,bound_map]=cell_seg(im);
                           
                           bw_dil = cast(bw_dil*255,'uint8');
                           imwrite(bw_dil,fsegname);
                           
                           if (saveoverlaid)
                               fovlname = fout_slim_overlaid(outdir,ff,tt,ii,chh,cc,rr,zz);
                               %Conver the image to the rgb as well
                               im = (im-min_phaseval)*255.0/(max_phaseval-min_phaseval);
                               %Chop out saturated values
                               im = im.*(im<255.0)+255.0*(im>=255.0);
                               im = im.*(im>=0);
                               %Convert to gray
                               rgbim = repmat(im,[1 1 3]);
                               rgbim = cast(rgbim,'uint8');

                               %Set all pixels on the boundaries to red
                               [nrows,ncols] = size(im);
                               figure(1)
                               imagesc(im);
                               hold on;
                               for cellidx = 1:ncells
                                   boundary = bound_map{cellidx};
                                   idx = (boundary(:,2)-1)*nrows+boundary(:,1);
                                   for channel = 1:3
                                        temp = rgbim(:,:,channel);
                                        if (channel==1)
                                            temp(idx)=255;
                                        else
                                            temp(idx)=0;
                                        end
                                        rgbim(:,:,channel) = temp;
                                   end
                                   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
                               end

                               imwrite(rgbim,fovlname);
                           end
                           save(fmatname,'ncells','bound_map');
                    end
               end
           end
      end
 end
 
