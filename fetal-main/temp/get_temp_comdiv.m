% UCSF
% 2022
% Jacob Ellison
% Calculating Relative Fetal Temperature
%
% Supported in matlab_2019a
% When running this script, you must have added it to the Matlab path and run script from the data directory
% You must change 'path_to_home_dir' at the end of the script to the location where you downloaded the 'temp' folder
%
%
% -------------------------------------------------------------------------------------------------------
t_file = dir('P*.7_*');


tempd_1=[];mag1=[];
tempd_2=[];mag2=[];
mask1=[];mask2=[];

% -- parsing p-files an calling com_div.m -- 
i=1
while isempty(tempd_2)
    t1_file = t_file(i).name;
    num1=(split(t_file(i).name, '_')); %need char to string for newer matlab version
    time1=str2num(num1{end});
    timepoint1=num1{end};
    pfile1=num1{1};
    scan1=split(pfile1,'P');
    scan1=str2num(scan1{end});
    for k=1:(length(t_file))
        num2=(split(t_file(k).name, '_'));
        time2=str2num(num2{end});
        timepoint2=num2{end};
        pfile2=num2{1};
        scan2=split(pfile2,'P');
        scan2=str2num(scan2{end});
        if (time1 == time2 | time1 ==(time2-1) | time1==time2-2) & (min(pfile1==pfile2)==0) & (scan1==scan2-512)

            sprintf(t_file(i).name)
            sprintf(t_file(k).name)
            [tempd,mag,m]=com_div(t_file(i).name,t_file(k).name);

            if isempty(tempd_1)
                scan1_tp = timepoint2;
                tempd_1=tempd;
                mag1=mag;
                mask1=m;
            else
                scan2_tp = timepoint1;
                tempd_2=tempd;
                mag2=mag;
                mask2=m;
            end

        end
    end
    i=i+1;    
end

%-- get scan duration -- 
dt2=datetime([2020, str2num(scan2_tp(1:2)),str2num(scan2_tp(3:4)),str2num(scan2_tp(5:6)), str2num(scan2_tp(7:8)), 0 ]);
dt1=datetime([2020, str2num(scan1_tp(1:2)),str2num(scan1_tp(3:4)),str2num(scan1_tp(5:6)), str2num(scan1_tp(7:8)), 0 ]);
duration = dt2-dt1

%-- get number of t2-weighted SSFSE images aquired -- 
scans = get_num_scans

% -- plotting -- 
figure;
subplot(121)
imagesc(tempd_1)
title('raw temp change at first scan')
subplot(122)
imagesc(tempd_2)
title('raw temp change at second scan')

% -- get final mask (overlap)--
mask_final = mask1.*mask2;

% -- Otsus thresholding on Movement (not used) -- 
Tot_dif = (tempd_2.*mask_final)-(tempd_1.*mask_final);
movement = mag2./(max(max(mag2))) - mag1./(max(max(mag1)));
thresh=ones(size(movement)); %set thresh to all ones in the case where this is not used

% -- sepereate gluteal ROI (small) fron brain ROI (big) --
mask_brain = bwareafilt(logical(mask_final),1, 'largest');
mask_back = bwareafilt(logical(mask_final),1, 'smallest');

% -- set zeros outside of ROIs to nans to remove later when averaging  -- 
t1_brain_temp = mask_brain.* thresh .* tempd_1;
t1_brain_temp(t1_brain_temp==0)=NaN;
t1_back_temp= mask_back .* thresh .* tempd_1;
t1_back_temp(t1_back_temp==0)=NaN;
t2_brain_temp = mask_brain .* thresh .* tempd_2;
t2_brain_temp(t2_brain_temp==0)=NaN;
t2_back_temp= mask_back .* thresh .* tempd_2;
t2_back_temp(t2_back_temp==0)=NaN;

% -- average phase in ROIs at t1 and t2 --
t1_mother_back_fat_avg_dtemp=mean(t1_back_temp(:),'omitnan')
t1_fetal_avg_dtemp=mean(t1_brain_temp(:),'omitnan')

t2_mother_back_fat_avg_dtemp=mean(t2_back_temp(:),'omitnan')
t2_fetal_avg_dtemp=mean(t2_brain_temp(:),'omitnan')

% -- calculate t1 and t2 relative shifts -- 
t2_shift = (t2_fetal_avg_dtemp-t2_mother_back_fat_avg_dtemp)
t1_shift = (t1_fetal_avg_dtemp-t1_mother_back_fat_avg_dtemp)

% -- calculate total relative phase --
z_phase=(t2_fetal_avg_dtemp-t2_mother_back_fat_avg_dtemp) - (t1_fetal_avg_dtemp-t1_mother_back_fat_avg_dtemp)


% -- plotting -- 
figure
subplot(3,4,[1 2])
imagesc(mag2);colormap gray;title('Magnitude Image'); ; axis off
subplot(3,4,[3 4])
imagesc(Tot_dif );colormap ;title('Raw Temperature Difference'); ; axis off
subplot(3,4,[5 6])
imagesc(movement);title('Movement Corruption'); ; axis off
subplot(3,4,[7 8])
imagesc(thresh);title('Movement Resulting in > 5% Signal Change'); ; axis off
subplot(3,4,[10 11])
imagesc(thresh.*Tot_dif);title('Movement Thresholded Tempurature Difference'); ; axis off

% -- alternate -- 
a =t2_mother_back_fat_avg_dtemp- t1_mother_back_fat_avg_dtemp;
figure;imagesc(phase2temp((t2_brain_temp-t1_brain_temp - a),12,3));title('Mapped relative dT'); colormap gray; axis off
b=(t2_brain_temp-t1_brain_temp - a);
mean(b(:), 'omitnan'),median(b(:), 'omitnan')


% -- save ROIs --
d=pwd;
s=split(d,'/');
cur = s{end};
new_fold=join(['path_to_home_dir/temp/records_temp/',cur, '/']);
mkdir(new_fold);
file_mask = join(['path_to_home_dir/temp/records_temp/',cur, '/',cur,'_mask.mat']);
save(file_mask,'mask_final');

f = labeloverlay(mag2,(mask_brain+ 2.*mask_back),'Colormap','jet','Transparency',0.7,'Includedlabels',[1,2]);
fig = figure;imagesc(f);title('Overlay T2');axis off; axis equal;
file_ol2 = join(['path_to_home_dir/temp/records_temp/',cur, '/',cur,'_t2.png']);
saveas(fig,file_ol2)

f = labeloverlay(mag1,(mask_brain+ 2.*mask_back),'Colormap','jet','Transparency',0.7,'Includedlabels',[1,2]);
file_ol1 = join(['path_to_home_dir/temp/records_temp/',cur, '/',cur,'_t1.png']);
fig = figure;imagesc(f);title('Overlay T1');axis off; axis equal;
saveas(fig,file_ol1)


%-- final metrics saved to csv --
subj = {join(['subj',cur])};
ROI = {cur};

duration = [minutes(duration)];
scans = [scans];
t1mom = [t1_mother_back_fat_avg_dtemp];
t1fetal	= [t1_fetal_avg_dtemp];
t2mom = [t2_mother_back_fat_avg_dtemp];
t2fetal	= [t2_fetal_avg_dtemp];
t2shift	= [t2_shift];
t1shift	= [t1_shift];
relativePhase = [z_phase];
relativeTemp = [relativePhase/0.0963];
% -- derivation for 0.0963 as conversion to temperature through PRF -- 
% B0 = 3
% TE = 12
% B = B0; % units [T]
% TE = TE*10^(-3); % units [s]
% gamma_bar = 42.575*10^(6); % gyromagnetic ratio units [Hz/T]
% gamma_bar_rad = gamma_bar * (2*pi); %units [rad/s*T]
% alpha = -0.01; % thermal coefficient units [ppm/degC]
% alpha = 0.01*10^(-6); % thermal coefficient units [/degC]
% denom = gamma_bar_rad*alpha; % units [rad/s*T*degC]
% denom = denom*B*TE; % units [rad*s*T/s*T*degC] => [rad / degC]
% -----------------------------------------------------------
x = input('How bad was movement between timepoints (scale 0-10)?');
movement = [x];

T = table(subj,movement,duration,scans,t1mom,t1fetal,t2mom,t2fetal,t2shift,t1shift,relativePhase,relativeTemp, 'RowNames',ROI);
file_csv = join(['path_to_home_dir/temp/records_temp/',cur, '/',cur,'_metrics.csv']);
writetable(T,file_csv);
close all