function putallimages (handle, images, slices, frames)
%PUTALLIMAGES  write whole or partial images to an open MINC file.
%
%  putallimages (handle, images)
%     write the matrix 'images' of size [nvox, nslices, nframes] to
%     the open MINC file specified by 'handle'. The size of the 
%     matrix must match the size of the open image volume.
%
%  putallimages (handle, images, slices, frames)
%     'slices' and 'frames' are vectors containing indices to
%     the slices and frames to be written. The size of the input
%     matrix must be [nvoxels, > max(slices), > max(frames)].
%     If 'slices' is a scalar, the input matrix can be of size
%     [nvoxels, > max(frames)].
%
% last modified $Date: 2003/06/26 02:51:18 $
% by            $Author: warnking $


% get some info about the data, if not provided: write entire matrix by default
if nargin < 3,
   nslices = size(images,2);
   nframes = size(images,3);
   slices = [1:nslices];
   frames = [1:nframes];
end;


% set the empirically deduced limit on the number of frames
max_frames = 150;

%check input
if ((length(size(images)) == 2) & (length(slices)==1)) % 'old' fashion of calling putallimages from inside a slices loop
   % pad image matrix to have data in the requested slice (quick fix)
   nvoxels = size(images,1);
   nframes = size(images,2);
   % create a slice dimension in image matrix
   images = reshape(images,[nvoxels, 1, nframes]); 
   % repeat data along the slice dimension so the new routine is happy
   images = images(:,ones(slices),:);
elseif (size(images,2) < max(slices))
   error('image matrix passed does not have enough slices');
end;
if size(images,3) < max(frames)
   error('image matrix passed does not have enough frames');
end;

% figure out how many times to call putimages
num_frames = size(frames,2);
num_calls = ceil(num_frames/max_frames);

for k = slices,
   
   % initialize start frame
   start_frame=1;
   
   % call putimages incrementally
   for i=1:num_calls
      
      end_frame = min(i*max_frames, num_frames);
      putimages(handle,images(:,k,start_frame:end_frame),k,frames(start_frame:end_frame));
      start_frame = end_frame+1;
      
   end;
   
end;
