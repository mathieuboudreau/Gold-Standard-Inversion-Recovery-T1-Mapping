function [ret, rname_out, iname_out] = convert_mp_ri_ym(mag_img,ph_img,varargin)
% [ret, rname_out, iname_out, pname_out] = convert_mp_ri(mbase,pbase,[prange,filter_phase])
%
% Convert the input magnitude/phase data to real/imaginary data
% Also output the actual phase images.
%
% IN: mbase - the file name base of the magnitude minc file
%     pbase - the file name base of the corresponding phase minc file
%
%     prange (optional) - range to be used for conversion of phase
%     filter_phase (optional) - flag for median filtering of phase map prior to Re/Im operation
%   The order of input is VERY important!
%
% OUT: ret (= 1 if successful)
%      rname_out (= mbase '_ll_t1_real.mnc') - the real minc file
%      iname_out (= pbase '_ll_t1_imag.mnc') - the imaginary minc file
%      pname_out (= pbase '_phase.mnc')      - the actual phase values minc file
%      rfname_out (= mbase '_ll_t1_real_filt.mnc') - the real minc file, after phase filtering
%      ifname_out (= pbase '_ll_t1_imag_filt.mnc') - the imaginary minc file, after phase filtering
%      pfname_out (= pbase '_phase_filt.mnc')      - the actual phase values minc file, after phase filtering
%
%	version 0.3

%------------------------------------------------------------------------
%
%    Coded by Bojana Stefanovic, c. 2002-2003?
%
%    Mod. July/06: allow optional specification of phase image range (Ives Levesque)
%    Mod. Oct/06:  force specification of entire mag/phase image names
%    Mod. Oct/07:  use proper flip angle value (20 degs)
%    Mod. Feb/11:  added flag for optional median filtering of phase data 
%    Mod. Feb/12: by Yuhan Ma: get rid of reading time frames.
%    getallimages.m is also modified accordingly
%
%------------------------------------------------------------------------

mbase = mag_img(1:findstr('.mnc',mag_img)-1);
pbase = ph_img(1:findstr('.mnc',ph_img)-1);

% set the range for rescaling phase values to radians on the range of [-pi, pi]
default_range = 4096;
if nargin > 2
    prange = varargin{1};
else
    prange = default_range;
end

if nargin > 3
    filter_phase = varargin{2};
    
else
    filter_phase = 0;
end

% set the time to first EX and the interval btw EXs
% HARD-CODED
start = 0.015+0.00079;
step = 0.495;
% set the corrected flip angle
%fa = 18.91;    % This is for the old Vision...  update?
fa = 20;	% updated for Sonata

% open images and get info
%magn = strcat(mbase,'.mnc');
m = openimage(mag_img,'r');
% nframes = getimageinfo(m,'NumFrames');

nframes = 1;
nslices = getimageinfo(m,'NumSlices');
height  = getimageinfo(m,'ImageHeight');
width   = getimageinfo(m,'ImageWidth');

%phase = strcat(pbase,'.mnc');
p = openimage(ph_img,'r');
% for data with time frames: use the following
% phase_meas  = getallimages_ym(p,1,1:nframes);
% otherwise, no time frames,use:
phase_meas  = getallimages_ym(p,1);

% MinMax returns the minimum and maximum value for the whole volume
minmax = getimageinfo(p,'MinMax');

% check dynamic range of phase image
img_range = abs(minmax(1) - minmax(2));
if img_range > prange
    warning(sprintf('Specified range value (%d) is not great enough to cover phase image dynamic range (%d). Errors may result', prange, img_range));
end


% intialize output data
real_data=NaN*ones(height*width,nframes);
imag_data=NaN*ones(height*width,nframes);
if filter_phase
    realf_data=NaN*ones(height*width,nframes);
    imagf_data=NaN*ones(height*width,nframes);
end


% write the real/imag data to new files
% rname_out = strcat(mbase,'_ll_t1_real.mnc');
% ym:


rname_out = strcat(mbase,'_real.mnc');

fprintf('Writing %s\n',rname_out);
%r_out = newimage(rname_out, [nframes nslices height width], mag_img);

r_out = newimage(rname_out, [nframes nslices height width], mag_img, [], [-32768 32767]);

% iname_out = strcat(pbase,'_ll_t1_imag.mnc');
iname_out = strcat(pbase,'_imag.mnc');
fprintf('Writing %s\n',iname_out);
%i_out = newimage(iname_out, [nframes nslices height width], mag_img);
i_out = newimage(iname_out, [nframes nslices height width], mag_img, [], [-32768 32767]);

pname_out = strcat(pbase,'_phase.mnc');
fprintf('Writing %s\n',pname_out);
%p_out = newimage(pname_out, [nframes nslices height width], mag_img);
p_out = newimage(pname_out, [nframes nslices height width], mag_img, [], [-32768 32767]);

if filter_phase
    % write the filtered real/imag data to new files
    rfname_out = strcat(mbase,'_ll_t1_real_filt.mnc');
    fprintf('Writing %s\n',rfname_out);
    %r_out = newimage(rname_out, [nframes nslices height width], mag_img);
    rf_out = newimage(rfname_out, [nframes nslices height width], mag_img, [], [-32768 32767]);

    ifname_out = strcat(pbase,'_ll_t1_imag_filt.mnc');
    fprintf('Writing %s\n',ifname_out);
    %i_out = newimage(iname_out, [nframes nslices height width], mag_img);
    if_out = newimage(ifname_out, [nframes nslices height width], mag_img, [], [-32768 32767]);

    pfname_out = strcat(pbase,'_phase_filt.mnc');
    fprintf('Writing %s\n',pfname_out);
    %p_out = newimage(pname_out, [nframes nslices height width], mag_img);
    pf_out = newimage(pfname_out, [nframes nslices height width], mag_img, [], [-32768 32767]);
end

% convert M/P to R/I for every slice
for slice=1:nslices
  
  % get data with frames
%   magn = getimages(m,slice,1:nframes);
%   phase_meas = getimages(p,slice,1:nframes);
  % ym: no frames
  magn = getimages(m,slice);
  phase_meas = getimages(p,slice);

  % calculate actual phase data (relative to the first frame)
  phase = (phase_meas - zeros(size(phase_meas)))*(2*pi/prange);

  % subtract pi from the phase data
  phase = phase - pi;

  % optional median filter the phase map for the slice
  if filter_phase
    phase_filt = zeros(size(phase));
    x_filt = 3;
    y_filt = 3;
    for f = 1:nframes  % loop over frames to filter each
      tmp = reshape(phase(:,f),width,height);
      tmp = medfilt2(tmp,[x_filt,y_filt]);
      phase_filt(:,f) = reshape(tmp,height*width,1);
    end
  end

  
  % calculate real and imaginary parts
  real_data = magn.*cos(phase);
  imag_data = magn.*sin(phase);      

  % write out the current frame of real/imaginary data
  putallimages(r_out,real_data,slice,1:nframes);
  putallimages(i_out,imag_data,slice,1:nframes);

  % write out the current frame of phase data
  putallimages(p_out,phase,slice,1:nframes);

  clear phase_meas phase real_data imag_data;

  if filter_phase
    % calculate real and imaginary parts
    realf_data = magn.*cos(phase_filt);
    imagf_data = magn.*sin(phase_filt);      

    % write out the current frame of real/imaginary data
    putallimages(rf_out,realf_data,slice,1:nframes);
    putallimages(if_out,imagf_data,slice,1:nframes);

    % write out the current frame of phase data
    putallimages(pf_out,phase_filt,slice,1:nframes);

    clear phase_filt realf_data imagf_data;
  end
  
  clear magn;
  
end

% write start (t to first EX) and step (t btw EXs) into the mincheaders
unix(['minc_modify_header -dinsert time:start=',num2str(start),' ',rname_out]);
unix(['minc_modify_header -dinsert time:start=',num2str(start),' ',iname_out]);

unix(['minc_modify_header -dinsert time:step=',num2str(step),' ',rname_out]);
unix(['minc_modify_header -dinsert time:step=',num2str(step),' ',iname_out]);

unix(['minc_modify_header -dinsert acquisition:delay_time=',num2str(step),' ',rname_out]);
unix(['minc_modify_header -dinsert acquisition:delay_time=',num2str(step),' ',iname_out]);

% write the corrected flip angle to the mincheader
unix(['minc_modify_header -dinsert acquisition:flip_angle=',num2str(fa),' ',rname_out]);
unix(['minc_modify_header -dinsert acquisition:flip_angle=',num2str(fa),' ',iname_out]);

% close all files
closeimage(m);
closeimage(p);
closeimage(r_out);
closeimage(i_out);
closeimage(p_out);


if filter_phase
    % write start (t to first EX) and step (t btw EXs) into the mincheaders
    unix(['minc_modify_header -dinsert time:start=',num2str(start),' ',rfname_out]);
    unix(['minc_modify_header -dinsert time:start=',num2str(start),' ',ifname_out]);

    unix(['minc_modify_header -dinsert time:step=',num2str(step),' ',rfname_out]);
    unix(['minc_modify_header -dinsert time:step=',num2str(step),' ',ifname_out]);

    unix(['minc_modify_header -dinsert acquisition:delay_time=',num2str(step),' ',rfname_out]);
    unix(['minc_modify_header -dinsert acquisition:delay_time=',num2str(step),' ',ifname_out]);

    % write the corrected flip angle to the mincheader
    unix(['minc_modify_header -dinsert acquisition:flip_angle=',num2str(fa),' ',rfname_out]);
    unix(['minc_modify_header -dinsert acquisition:flip_angle=',num2str(fa),' ',ifname_out]);

    closeimage(rf_out);
    closeimage(if_out);
    closeimage(pf_out);
    
end

ret = 1;







