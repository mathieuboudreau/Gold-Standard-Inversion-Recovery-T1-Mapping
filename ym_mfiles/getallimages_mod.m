function images = getallimages (handle, slices, frames)
%GETALLIMAGES  Retrieve whole or partial images from a MINC file.
%
%  images = getallimages (fname)
%     read the entire file specified by fname
%
%  images = getallimages (fname, slices, frames)
%     'slices' and 'frames' are vectors containing indices to
%     the slices and frames to be read. The size of the return
%     matrix is [nvoxels, max(slices), max(frames)].
%     If 'slices' is a scalar, the returned matrix is of size
%     [nvoxels, max(frames)].
%
%  images = getallimages (handle, ...)
%     pass in the handle of a file that is already open, instead
%     of a filename.
%
% last modified $Date: 2003/06/26 02:51:18 $
% by            $Author: warnking $

% if we got passed a filename, not a handle
% if ischar(handle), 
%    handle = openimage(handle,'r');
%    close_image = 1;
% else
%    close_image = 0;
% end;

% get some info about the file, if not provided: read all files by default
if nargin == 1,
   nslices = getimageinfo(handle,'NumSlices');
   
   % nframes = getimageinfo(handle,'NumFrames');
   nframes = 1;
   slices = [1:nslices];
   % frames = [1:nframes];
end;

% set the empirically deduced limit on the number of frames
max_frames = 150;

% figure out how many times to call getimages

% num_frames = size(frames,2);

% ym: no time frames
num_frames = 1;
num_calls = ceil(num_frames/max_frames)

for k = slices
   
   % initialize start frame
   start_frame=1;
   
%    % call getimages incrementally
%    for i=1:num_calls
%       
%       end_frame = min(i*max_frames, num_frames);
%       images(:,k,start_frame:end_frame) = getimages(handle, k, frames(start_frame:end_frame));
%       size(images)
%       start_frame = end_frame+1;
%       
%    end;
   
end;

% check for 'old' fashion of calling getallimages from inside a slices loop
if (length(slices)==1) 
   images = squeeze(images(:,k,:));
end;

% close the image, if we opened it
if close_image,
   closeimage(handle);
end;