%supported in matlab_2019a

t_file = dir('P*.7_*');


tempd_1=[];mag1=[];
tempd_2=[];mag2=[];
mask1=[];mask2=[];
% for i =1:(length(t_file))
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
            % [tempd,mag,m]=com_div(t_file(i).name,t_file(k).name);
            tempd=1;
            mag=1;
            m=1;



            if isempty(tempd_1)
                scan1_tp = timepoint2;
                tempd_1=tempd;
                % mag1=mag;
                % mask1=m;
            else
                scan2_tp = timepoint1;
                tempd_2=tempd;
                % mag2=mag;
                % mask2=m;
            end

        end
    end
    i=i+1;    
end

dt2=datetime([2020, str2num(scan2_tp(1:2)),str2num(scan2_tp(3:4)),str2num(scan2_tp(5:6)), str2num(scan2_tp(7:8)), 0 ]);
dt1=datetime([2020, str2num(scan1_tp(1:2)),str2num(scan1_tp(3:4)),str2num(scan1_tp(5:6)), str2num(scan1_tp(7:8)), 0 ]);

duration = dt2-dt1
scan1_tp
scan2_tp
% num2
% % print('num scans:')
% scans = get_num_scans
% % tempd_1 = tempd_1 .* -1;
% figure;
% subplot(121)
% imagesc(tempd_1)
% title('raw temp change at first scan')
% subplot(122)
% imagesc(tempd_2)
% title('raw temp change at second scan')

% mask_final = mask1.*mask2;

% Tot_dif = (tempd_2.*mask_final)-(tempd_1.*mask_final);
% movement = mag2./(max(max(mag2))) - mag1./(max(max(mag1)));

% thresh=ones(size(movement));
% % thresh(movement>0.05)=1; % determine why this range
% % thresh(movement<-0.05)=1;
% % tempdiff=phase2temp((rpdif_cmb.*thresh),12,3);
% % tempdiff=phase2temp((up.*thresh),12,3);


% % back=bp_2-bp_1;
% % tempdiff_back=phase2temp((back.*thresh),12,3); % determine the artifacts and what could be the causes of hotspots

% % fff=tempdiff>0;
% % td_brain=mean(mean(mean(tempdiff(fff))))
% % fetal_brain_avg_dtemp=mean(tempdiff(:),'omitnan')


% % fff=tempdiff_back>0;
% % td_back=mean(mean(mean(tempdiff_back(fff))))
% % mother_back_fat_avg_dtemp=mean(tempdiff_back(:),'omitnan')

% % %==========final mask ==========
% % figure;
% % % subplot(121);imagesc(movement);colorbar
% % imagesc(Tot_dif)
% % h = imfreehand();
% % mask_brain = h.createMask();
% % mask_brain=double(mask_brain);
% % mask_brain(mask_brain==0)=NaN;
% % close
% % % T = graythresh(mag2.*mask_brain);
% % % mask_brain = (mask_brain.*(mag2>T));

% % figure;
% % % subplot(121);imagesc(movement);colorbar
% % imagesc(Tot_dif)
% % h = imfreehand();
% % mask_back = h.createMask();
% % mask_back=double(mask_back);
% % mask_back(mask_back==0)=NaN;
% % close
% % % brain_temp= phase2temp(( mask_brain .* thresh .* Tot_dif),12,3);
% % % back_temp= phase2temp((mask_back .* thresh .* Tot_dif),12,3);
% % %========final mask ===========

% mask_brain = bwareafilt(logical(mask_final),1, 'largest');
% mask_back = bwareafilt(logical(mask_final),1, 'smallest');
% % CC = bwconncomp(mask_final);
% % final_brain_mask=mask_final;
% % final_back_mask=mask_final;
% % numOfPixels = cellfun(@numel,CC.PixelIdxList);
% % [unused,indexOfMax] = max(numOfPixels);
% % final_brain_mask(CC.PixelIdxList{indexOfMax}) = 0;
% % figure;imagesc(final_brain_mask);



% t1_brain_temp = mask_brain.* thresh .* tempd_1;
% t1_brain_temp(t1_brain_temp==0)=NaN;
% t1_back_temp= mask_back .* thresh .* tempd_1;
% t1_back_temp(t1_back_temp==0)=NaN;
% t2_brain_temp = mask_brain .* thresh .* tempd_2;
% t2_brain_temp(t2_brain_temp==0)=NaN;
% t2_back_temp= mask_back .* thresh .* tempd_2;
% t2_back_temp(t2_back_temp==0)=NaN;



% t1_mother_back_fat_avg_dtemp=mean(t1_back_temp(:),'omitnan')
% t1_fetal_avg_dtemp=mean(t1_brain_temp(:),'omitnan')


% t2_mother_back_fat_avg_dtemp=mean(t2_back_temp(:),'omitnan')
% t2_fetal_avg_dtemp=mean(t2_brain_temp(:),'omitnan')


% % mother_back_fat_avg_dtemp=mean(back_temp(:),'omitnan');
% % fetal_avg_dtemp=mean(brain_temp(:),'omitnan');

% t2_shift = (t2_fetal_avg_dtemp-t2_mother_back_fat_avg_dtemp)
% t1_shift = (t1_fetal_avg_dtemp-t1_mother_back_fat_avg_dtemp)
% % z_phase_im=(t2_brain_temp-t2_back_temp) - (t1_brain_temp-t1_back_temp);
% z_phase=(t2_fetal_avg_dtemp-t2_mother_back_fat_avg_dtemp) - (t1_fetal_avg_dtemp-t1_mother_back_fat_avg_dtemp)
% % rel_tem =z_phase/denom


% % figure
% % subplot(121)
% % imagesc(brain_temp)
% % title('Total fetal brain temp change')
% % subplot(122)
% % imagesc(back_temp)
% % title('Total mother back fat temp change')



% figure
% subplot(3,4,[1 2])
% % str = ({'fetal brain avg temp diff', strcat(num2str(fetal_brain_avg_dtemp), ' ^{\circ} C')})
% imagesc(mag2);colormap gray;title('Magnitude Image'); ; axis off
% subplot(3,4,[3 4])
% % str = ({'fetal brain avg temp diff', strcat(num2str(fetal_brain_avg_dtemp), ' ^{\circ} C')})
% imagesc(Tot_dif );colormap ;title('Raw Temperature Difference'); ; axis off

% subplot(3,4,[5 6])
% % str = ({'fat pad avg temp diff', strcat(num2str(mother_back_fat_avg_dtemp), ' ^{\circ} C')})
% imagesc(movement);title('Movement Corruption'); ; axis off

% subplot(3,4,[7 8])
% % str = ({'fat pad avg temp diff', strcat(num2str(mother_back_fat_avg_dtemp), ' ^{\circ} C')})
% imagesc(thresh);title('Movement Resulting in > 5% Signal Change'); ; axis off

% subplot(3,4,[10 11])
% % str = ({'fat pad avg temp diff', strcat(num2str(mother_back_fat_avg_dtemp), ' ^{\circ} C')})
% imagesc(thresh.*Tot_dif);title('Movement Thresholded Tempurature Difference'); ; axis off


% % figure
% % subplot(221)
% % imagesc(mag1);colormap gray;title('magnitude GRE t1'); axis off; colorbar 
% % subplot(222)
% % imagesc(mag2);colormap gray;title('magnitude GRE t2'); axis off; colorbar

% % subplot(223)
% % imagesc(rmdif);colormap gray;title('total mag diff (normalized)'); axis off; colorbar
% % subplot(224)
% % imagesc(thresh);title('movement threshold'); axis off; colorbar

% % figure
% % subplot(221)
% % imagesc(pdif1_cmb);colormap gray;title('hpf phase diff GRE t1 (rad)'); axis off; colorbar
% % subplot(222)
% % imagesc(pdif2_cmb);colormap gray;title('hpf phase diff GRE t2 (rad)'); axis off; colorbar
% % subplot(223)
% % imagesc(un_phase1);colormap gray;title('phase diff GRE t1 (rad)'); axis off; colorbar
% % subplot(224)
% % imagesc(un_phase2);colormap gray;title('phase diff GRE t2 (rad)'); axis off; colorbar

% % %dicomwrite(mag1,'/working/lupolab/jacob/fetal/temp/segment_test/brain.dcm')

% % %/working/lupolab/jacob/fetal/temp/segment_test/brain_seg.nii.gz

% a =t2_mother_back_fat_avg_dtemp- t1_mother_back_fat_avg_dtemp;
% figure;imagesc(phase2temp((t2_brain_temp-t1_brain_temp - a),12,3));title('Mapped relative dT'); colormap gray; axis off
% b=(t2_brain_temp-t1_brain_temp - a);
% mean(b(:), 'omitnan'),median(b(:), 'omitnan')

% d=pwd;
% s=split(d,'/');
% cur = s{end};
% new_fold=join(['/working/lupolab/jacob/fetal/temp/records_temp/',cur, '/']);
% mkdir(new_fold);
% file_mask = join(['/working/lupolab/jacob/fetal/temp/records_temp/',cur, '/',cur,'_mask.mat']);
% save(file_mask,'mask_final');
% % lbl_fin = (mask_brain+ 2.*mask_back);
% % lbl_fin(lbl_fin=1)='brain';
% f = labeloverlay(mag2,(mask_brain+ 2.*mask_back),'Colormap','jet','Transparency',0.7,'Includedlabels',[1,2]);
% fig = figure;imagesc(f);title('Overlay T2');axis off; axis equal;
% file_ol2 = join(['/working/lupolab/jacob/fetal/temp/records_temp/',cur, '/',cur,'_t2.png']);
% saveas(fig,file_ol2)

% f = labeloverlay(mag1,(mask_brain+ 2.*mask_back),'Colormap','jet','Transparency',0.7,'Includedlabels',[1,2]);
% file_ol1 = join(['/working/lupolab/jacob/fetal/temp/records_temp/',cur, '/',cur,'_t1.png']);
% fig = figure;imagesc(f);title('Overlay T1');axis off; axis equal;
% saveas(fig,file_ol1)

% subj = {join(['subj',cur])};
% ROI = {cur};

% duration = [minutes(duration)];
% scans = [scans];
% t1mom = [t1_mother_back_fat_avg_dtemp];
% t1fetal	= [t1_fetal_avg_dtemp];
% t2mom = [t2_mother_back_fat_avg_dtemp];
% t2fetal	= [t2_fetal_avg_dtemp];
% t2shift	= [t2_shift];
% t1shift	= [t1_shift];
% relativePhase = [z_phase];
% relativeTemp = [relativePhase/0.0963];
% x = input('How bad was movement between timepoints (scale 0-10)?');
% movement = [x];

% T = table(subj,movement,duration,scans,t1mom,t1fetal,t2mom,t2fetal,t2shift,t1shift,relativePhase,relativeTemp, 'RowNames',ROI);
% file_csv = join(['/working/lupolab/jacob/fetal/temp/records_temp/',cur, '/',cur,'_metrics.csv']);
% writetable(T,file_csv);
% close all