
function [t1,mag_cmb, mask_fin] = com_div(pfile,pfile1);
    % This function takes two p-files and returns the phase evolved between them for fetal brain and mother Gluteal tissue ROIs
    % ---- recon ----
    [k hdr] = read_MR_rawdata(pfile);
    hdr.data_acq_tab
    hdr.exam
    hdr.series

    size_p=size(k);
    ncoil=size_p(length(size_p));

    
    im = fftshift(ifft2(ifftshift(squeeze(k),1)),2);
    magni = abs(im);   % magnitude image
    phase = angle(im); % phase image
    magni = magni ./max(max(max(magni)));

    [k1 hdr1] = read_MR_rawdata(pfile1);
    size_p1=size(k1);
    ncoil1=size_p1(length(size_p1));

    im1 = fftshift(ifft2(ifftshift(squeeze(k1),1)),2);
    magni1 = abs(im1);   % magnitude image
    phase1 = angle(im1); % phase image
    magni1 = magni1 ./max(max(max(magni1)));

    % -- complex division --
    pase=(im1./im);

    % -- coil combination --
    mag_cmb = [zeros(size_p(1:2))];
    mag_cmb1 = [zeros(size_p(1:2))];
    cp_SUM = [zeros(size_p(1:2))];
    for i = 1:ncoil
        mag_cmb = mag_cmb + squeeze(magni(:,:,i)).^2;
        mag_cmb1 = mag_cmb1 + squeeze(magni1(:,:,i)).^2;
        cp_SUM = cp_SUM + pase(:,:,i);
    end
    pdc_coil1=angle(cp_SUM);
    mag_cmb=sqrt(mag_cmb);
    mag_cmb1=sqrt(mag_cmb1);

    % -- masking (hand drawn ROIs) -- 
    figure; colormap gray
    imagesc(mag_cmb1,[0,0.8]);
    h = imfreehand();
    mask_brain1 = h.createMask();
    mask_brain1=double(mask_brain1);
    close

    figure; colormap gray
    imagesc(mag_cmb1,[0,0.8]);
    h = imfreehand();
    mask_back1 = h.createMask();
    mask_back1=double(mask_back1);
    close
    mask1 = mask_brain1 + mask_back1;

    figure; colormap gray
    imagesc(mag_cmb,[0,0.8]);colormap gray
    h = imfreehand(); colormap gray
    mask_brain = h.createMask();
    mask_brain=double(mask_brain);
    close


    figure; colormap gray
    imagesc(mag_cmb,[0,0.8]); colormap gray
    h = imfreehand(); colormap gray
    mask_back = h.createMask();
    mask_back=double(mask_back);
    close
    mask = mask_brain + mask_back;

    %-- overlap in ROIs -- 
    mask_fin = mask .* mask1;

    %-- send abs of back final ROI -- 
    t1=abs(pdc_coil1.*mask_fin);
