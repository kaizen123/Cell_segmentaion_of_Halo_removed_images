
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
%    mkdir(matlab_dir);
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
f_dm_over_t=@(odir,t) sprintf('%s/total_drymass_at_%d.mat',odir,t); % FOV, TIME, Channel, Frame Number, PAT

%p=gcp;
%delete(p);
%p=parpool(8);
drymass = zeros(0,3); %Each row is a cell, the next is the total mass of slim, thresholded SLIM and the halo-removed SLIM
nbins = 200;
pixelratio = 3.2;%Pixel per micron
combined_slim_thr_hist = zeros(nbins,length(t));
combined_hr_hist = zeros(nbins,length(t));
update_dry_mass = 1;
override  = 0;
max_slim_thr_mass = 300;
 
if (update_dry_mass)
     for ff=f
          for tt=t  
              matfilename = f_dm_over_t(pwd,tt);
              drymass_over_time = zeros(0,3);
              if (~exist(matfilename)|override)
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
                                   minpixelnum = 10000; %Number of pixels in a cells
                                   maxpixelnum = 70000;
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
                                            if ((hr_dm<1000)&(length(curpixidxlist)>minpixelnum)&(length(curpixidxlist)<maxpixelnum))%If the cell is too small
                                                drymass(end+1,:)=[slim_dm slim_thres_dm hr_dm];
                                                drymass_over_time(end+1,:)=drymass(end,:);
                                            
                                            end
                                       end

                                  
                            end
                       end
                  end
              else
                  load(matfilename);
              end
              figure(1);
              plot(drymass_over_time(:,2),drymass_over_time(:,3),'+r');
              xlabel('Thresholded');
              ylabel('Halo removed');
              %drawnow;
              xhist = linspace(25,1000,nbins);
              xhist_slim_thr = linspace(25,max_slim_thr_mass,nbins);
             
                               
              cur_hr_total_dm_hist = hist(drymass_over_time(:,3),xhist);
              cur_slim_thr_dm_hist = hist(drymass_over_time(:,2),xhist_slim_thr);
              combined_slim_thr_hist(:,tt+1)=cur_slim_thr_dm_hist(:);
              combined_hr_hist(:,tt+1)=cur_hr_total_dm_hist(:);
              norm_combined_slim_thr_hist = combined_slim_thr_hist(2:end-1,:)./repmat(sum(combined_slim_thr_hist(2:end-1,:),1)+1e-6,[size(combined_slim_thr_hist(2:end-1,:),1) 1]);
              norm_combined_hr_hist = combined_hr_hist(2:end-1,:)./repmat(sum(combined_hr_hist(2:end-1,:),1)+1e-6,[size(combined_hr_hist(2:end-1,:),1) 1]);

              figure(2);
              subplot(221);
              imagesc(t,xhist_slim_thr(2:end-1),combined_slim_thr_hist(2:end-1,:));colormap jet;colorbar;title('Thresholded histogram');
              subplot(222);
              imagesc(t,xhist(2:end-1),combined_hr_hist(2:end-1,:));colormap jet;colorbar;title('HR histogram');
              subplot(223);
              imagesc(t,xhist_slim_thr(2:end-1),norm_combined_slim_thr_hist);colormap jet;colorbar;title('Norm. Thresholded histogram');
              subplot(224);
              imagesc(t,xhist(2:end-1),norm_combined_hr_hist);colormap jet;colorbar;title('Norm. HR histogram');
              %drawnow;

                  %save(matfilename,'drymass_over_time','cur_slim_thr_dm_hist','cur_hr_total_dm_hist',...
                  %    'norm_combined_slim_thr_hist','norm_combined_hr_hist');
               
                  
              %Compute the mean and the standard deviation of the drymass
              xhist_mat = repmat(xhist(:),[1 length(t)]);
              xhist_thr_mat = repmat(xhist_slim_thr(:),[1 length(t)]);
              mat1 = xhist_mat(2:end-1,:).*norm_combined_hr_hist;
              mat2 = xhist_thr_mat(2:end-1,:).*norm_combined_slim_thr_hist;
              mean_thr_dm = sum(mat2,1);
              mean_dm =sum(mat1,1);
              std_dm = sqrt(sum((xhist_mat(2:end-1,:) - repmat(mean_dm,[size(mat1,1) 1])).^2.*norm_combined_hr_hist,1));
              std_thr_dm = sqrt(sum((xhist_thr_mat(2:end-1,:) - repmat(mean_thr_dm,[size(mat1,1) 1])).^2.*norm_combined_slim_thr_hist,1));
              
              %std_thr_dm = std(mat1,1);
              %std_dm = std(mat2,1);
              
              
              figure(3);
              subplot(221);plot(t,mean_thr_dm);title('Mean Thres. SLIM ofver time');
              subplot(222);plot(t,mean_dm);title('Mean SLIM HR over time');
              subplot(223);plot(t,std_thr_dm);title('Std. Thres SLIM over time');
              subplot(224);plot(t,std_dm);title('Std. SLIM HR over time');
              
              
              end
          end

    
     save('All_dry_mass.mat','drymass','combined_slim_thr_hist','combined_hr_hist','norm_combined_slim_thr_hist','norm_combined_hr_hist',...
         'xhist','xhist_slim_thr');
else
    load('All_dry_mass.mat');%Load all dry mass values and also the histogram of drymass over the time
    %Display the scatter plot for all the mass
    %Compute the best fit coefficients
    slim_thr_mass = drymass(1:20:end,2);%A single cells is expected to appear 15 times...Lazy man approach!
    hr_mass = drymass(1:20:end,3);
    good_sample = find(slim_thr_mass<max_slim_thr_mass);
    slim_thr_mass=slim_thr_mass(good_sample);
    hr_mass = hr_mass(good_sample);
    
    alpha = sum(hr_mass)/sum(slim_thr_mass);
    disp(['Fitting coeff: ' num2str(alpha)]);
    %Draw all the data point
    figure(1);
    plot(slim_thr_mass,hr_mass,'.b'); 
    xlabel('Thresholded SLIM total drymass');
    ylabel('Halo removed total drymass');
    xmin = 0;
    xmax = max_slim_thr_mass;
    ymin = 0;
    ymax = 1020;
    axis([xmin xmax ymin ymax]);
    x_arr = linspace(xmin,xmax,100);
    f = @(x)(alpha*x);
    y_arr = f(x_arr);
    figure(1);hold on;
    plot(x_arr,y_arr,'-r','linewidth',3);
    %Compute R^2
    hr_mass_pred = f(slim_thr_mass);
    R2 = 1-norm(hr_mass-hr_mass_pred,'fro')^2/norm(hr_mass-mean(hr_mass),'fro')^2;
    disp(['R^2 values: ' num2str(R2)]);
    %Compute the minimum and maximum slope for the relation
    slope_arr = hr_mass./slim_thr_mass;
    min_slope = prctile(slope_arr,0.3);
    max_slope = prctile(slope_arr,99.7);
    fmin = @(x)(min_slope*x);
    fmax = @(x)(max_slope*x);
    disp(['Min slope: ' num2str(min_slope) ', Max slope: ' num2str(max_slope)]);
    figure(1);
    plot(x_arr,fmin(x_arr),'--r','linewidth',3);
    hold on
    plot(x_arr,fmax(x_arr),'--r','linewidth',3);
    
    %Display the histogram for drymass over the time. Make sure that we
    %reject the first and the last histogram band to ignore the extreme
    %values...
    
    time_interval = 22;
    figure(2);
    subplot(221);
    imagesc(t*time_interval,xhist_slim_thr(2:end-1),combined_slim_thr_hist(2:end-1,:));colormap jet;colorbar;title('Thresholded histogram');
    subplot(222);
    imagesc(t*time_interval,xhist(2:end-1),combined_hr_hist(2:end-1,:));colormap jet;colorbar;title('HR histogram');
    subplot(223);
    imagesc(t*time_interval,xhist_slim_thr(2:end-1),norm_combined_slim_thr_hist(2:end-1,:));colormap jet;colorbar;title('Norm. Thresholded histogram');
    subplot(224);
    imagesc(t*time_interval,xhist(2:end-1),norm_combined_hr_hist(2:end-1,:));colormap jet;colorbar;title('Norm. HR histogram');
    drawnow;

    
    
    
    
end

 %delete(p);
 
