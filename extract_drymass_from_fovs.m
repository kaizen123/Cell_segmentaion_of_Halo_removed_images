
%This code will extract all the drymasses of all the cells in the data and
%save it into a .mat file
%Author: Tan H. Nguyen
%University of Illinois at Urbana-Champaign
clc;
clear all;
close all;
outdir = 'D:\Hela_cell_time_laps_Feb_16th_2016\';
matlab_dir = strcat(outdir,'\mat_files\');
if (~exist(matlab_dir))
    mkdir(matlab_dir);
end
f=0:0;
t=[0:88];
chh=0;
ii=0;
r=0:24;
z=0:0;
c =0:24;
%File name for the Halo-removed image
fout_slim=@(odir,f,t,i,ch,c,r,z) sprintf('%s\\f%d_t%d_i%d_ch%d_c%d_r%d_z%d_SLIM.tif',odir,f,t,i,ch,c,r,z); % FOV, TIME, Channel, Frame Number, PAT
fout_slim_hr_ns=@(odir,f,t,i,ch,c,r,z) sprintf('%s\\f%d_t%d_i%d_ch%d_c%d_r%d_z%d_HR_NS.tif',odir,f,t,i,ch,c,r,z); % FOV, TIME, Channel, Frame Number, PAT

%File name for the segmentation image
fout_slim_seg=@(odir,f,t,i,ch,c,r,z) sprintf('%s\\f%d_t%d_i%d_ch%d_c%d_r%d_z%d_SEG.tif',odir,f,t,i,ch,c,r,z); % FOV, TIME, Channel, Frame Number, PAT
fout_slim_overlaid=@(odir,f,t,i,ch,c,r,z) sprintf('%s\\f%d_t%d_i%d_ch%d_c%d_r%d_z%d_OVL.tif',odir,f,t,i,ch,c,r,z); % FOV, TIME, Channel, Frame Number, PAT
fout_bound_mat=@(odir,f,t,i,ch,c,r,z) sprintf('%s\\f%d_t%d_i%d_ch%d_c%d_r%d_z%d_bound.mat',odir,f,t,i,ch,c,r,z); % FOV, TIME, Channel, Frame Number, PAT
f_dm_over_t=@(odir,t) sprintf('%s\\total_drymass_at_%d.mat',odir,t); % FOV, TIME, Channel, Frame Number, PAT

%p=gcp;
%delete(p);
%p=parpool(8);
drymass = zeros(0,3); %Each row is a cell, the next is the total mass of slim, thresholded SLIM and the halo-removed SLIM
nbins = 200;
pixelratio = 3.2;%Pixel per micron
combined_slim_thr_hist = zeros(nbins,length(t));
combined_hr_hist = zeros(nbins,length(t));
 for ff=f
      for tt=t  
          drymass_over_time = zeros(0,3);
          for cc=c 
               for rr = r
                    for zz = z
                            saveoverlaid = 0;
                            min_phaseval = -0.3;
                            max_phaseval = 2.5;

                           disp(['Processing r: ' num2str(rr) ', c: ' num2str(cc) ', t: ' num2str(tt)]);
                           fnsname = fout_slim_hr_ns(outdir,ff,tt,ii,chh,cc,rr,zz);
                           fslimname = fout_slim(outdir,ff,tt,ii,chh,cc,rr,zz); %Generate the filename for SLIM images
                           fsegname = fout_slim_seg(outdir,ff,tt,ii,chh,cc,rr,zz);
                           fmatname = fout_bound_mat(matlab_dir,ff,tt,ii,chh,cc,rr,zz); %Name of the mat file to be saved
                           bw_dil = imread(fsegname); %Read the bw segmented image
                           S = regionprops(im2bw(bw_dil),'PixelIdxList');
                           slim_im = imread(fslimname);
                           ncells = size(S,1); %Get the number of cells
                           slimim = single(imread(fslimname));
                           hrnsim = single(imread(fnsname));
                           for cellidx=1:ncells
                                %Get the current indices of all the pixels
                                curpixidxlist = S(cellidx).PixelIdxList;
                                %Total non-negative phase of slim
                                slim_total_phase = sum(slimim(curpixidxlist));
                                slim_thres_total_phase = sum(slimim(curpixidxlist).*(slimim(curpixidxlist)>0));
                                hr_total_phase = sum(hrnsim(curpixidxlist));
                                hr_dm = hr_total_phase*0.4235/pixelratio^2;
                                slim_dm = slim_total_phase*0.4235/pixelratio^2;
                                slim_thres_dm = slim_thres_total_phase*0.4235/pixelratio^2;
                                if ((hr_dm>50)&(hr_dm<1000))%If the cell is too small
                                    drymass(end+1,:)=[slim_dm slim_thres_dm hr_dm];
                                    drymass_over_time(end+1,:)=drymass(end,:);
                                end
                           end
                    end
               end
          end
          figure(1);
          plot(drymass_over_time(:,2),drymass_over_time(:,3),'+r');
          xlabel('Thresholded');
          ylabel('Halo removed');
          drawnow;
          xhist = linspace(50,1000,nbins);
          xhist_slim_thr = linspace(25,500,nbins);
          cur_hr_total_dm_hist = hist(drymass_over_time(:,3),xhist);
          cur_slim_thr_dm_hist = hist(drymass_over_time(:,2),xhist_slim_thr);
          matfilename = f_dm_over_t(pwd,tt);
          combined_slim_thr_hist(:,tt+1)=cur_slim_thr_dm_hist(:);
          combined_hr_hist(:,tt+1)=cur_slim_thr_dm_hist(:);
          norm_combined_slim_thr_hist = combined_slim_thr_hist./repmat(sum(combined_slim_thr_hist,1)+1e-6,[size(combined_slim_thr_hist,1) 1]);
          norm_combined_hr_hist = combined_hr_hist./repmat(sum(combined_hr_hist,1)+1e-6,[size(combined_hr_hist,1) 1]);
          
          figure(2);
          subplot(221);
          imagesc(t,xhist,combined_slim_thr_hist);colormap jet;colorbar;title('Thresholded histogram');
          subplot(222);
          imagesc(t,xhist_slim_thr,combined_hr_hist);colormap jet;colorbar;title('HR histogram');
          subplot(223);
          imagesc(norm_combined_slim_thr_hist);colormap jet;colorbar;title('Norm. Thresholded histogram');
          subplot(224);
          imagesc(norm_combined_hr_hist);colormap jet;colorbar;title('Norm. HR histogram');
          drawnow;
          
          save(matfilename,'drymass_over_time','cur_slim_thr_dm_hist','cur_hr_total_dm_hist',...
              'norm_combined_slim_thr_hist','norm_combined_hr_hist');
      end
      
 end
 save('All_dry_mass.mat','drymass','combined_slim_thr_hist','combined_hr_hist');
 %delete(p);
 
