% Load minc images, save in matrices
% Calculate real and imaginary parts from magnitude and phase images 
% Save real and imaginary parts in matrices

function complexConver(SUBJECT_ID,T1path, indeces_complex)

loadpath = T1path;

% Where to save the .mat 
%savename = [T1path 'data/' 'TestSingleSlice'];

savename = fullfile (T1path, 'data.mat');

currpath = pwd;
cd(loadpath)

%d: raw data
d=mincLoadAllSeries('.');
% indeces = [2,4,6,8];
% disp('plot d(5).imData')
% figure;
% imshow(d(5).imData,[]);

% indeces_complex
d=d(indeces_complex);

nbim = size(d,2)

% ym:

disp('plot d(1).imData')
figure;
imshow(d(1).imData,[]);
hold;


nbrow = size(d(1).imData, 1);
nbcol = size(d(1).imData, 2);
nbslices = size(d(1).imData, 3);
nbseries = length(d)/2;

% magDataTmp = zeros(nbrow,nbcol);
% phasDataTmp = zeros(nbrow,nbcol);
% realDataTmp = zeros(nbrow,nbcol);
% imagDataTmp = zeros(nbrow,nbcol);
% compexDataTmp = zeros(nbrow, nbcol);

magData = zeros(nbcol, nbrow, nbslices, nbseries);
phasData = zeros(nbcol, nbrow, nbslices, nbseries);
realData = zeros(nbcol, nbrow, nbslices, nbseries);
imagData = zeros(nbcol, nbrow, nbslices, nbseries);
data = zeros(nbcol, nbrow, nbslices, nbseries); %complex data

extra.tVec = zeros(1,nbseries); % One series corresponds to one TI



s=1;
kk = [1,3,5,7]

% from getDataMag:
% for k = 1:nbseries
% 	dataTmp = d(k).imData;
%    % sizedataTmp = size(dataTmp)
% 	dataTmp = double(squeeze(dataTmp));	
%    % sizeSqueezeDataTmp = size(dataTmp)
% 	for ss = 1:nbslice
%         % ym:
%         % sizedata = size(data)
% 		%data(:,:,ss,k) = dataTmp(:,:,3+(ss-1)*4)+i*dataTmp(:,:,4+(ss-1)*4); %SE
% 		%data(:,:,ss,k) = dataTmp(:,:,1+(ss-1)*4); %GRE
% 		data(:,:,ss,k) = reshape(dataTmp(:,:,ss), nbcol, nbrow); 
%         sizedatareshape = size(data)
% 	end
% 	extra.tVec(k) = d(k).inversionTime;
% end 


    for k=1:nbseries
        
            magDataTmp = d(s).imData;
            magDataTmp = double(squeeze(magDataTmp));
            phasDataTmp = d(s+1).imData;
            phasDataTmp = double(squeeze(phasDataTmp));
            
            for ss= 1:nbslices
                
                magData(:,:,ss,k)=reshape(magDataTmp(:,:,ss),nbcol, nbrow);
                phasData(:,:,ss,k)=reshape(phasDataTmp(:,:,ss),nbcol, nbrow);

                for m = 1:nbcol
                    for n= 1:nbrow
                        
                        realData(m,n,ss,k)=magData(m,n,ss,k)*cos(phasData(m,n,ss,k)*pi/180);
                        imagData(m,n,ss,k)=sqrt(-1)*magData(m,n,ss,k)*sin(phasData(m,n,ss,k)*pi/180);
%                         realDataTmp(m,n)=magDataTmp(m,n)*cos(phasDataTmp(m,n)*pi/180);
%                         imagDataTmp(m,n)=sqrt(-1)*magDataTmp(m,n)*sin(phasDataTmp(m,n)*pi/180);
%                         complexDataTmp(m,n)=realDataTmp(m,n)+imagDataTmp(m,n);

                          


                    end

                end
                
                data(:,:,ss,k)=realData(:,:,ss,k)+imagData(:,:,ss,k);
            end
            
           
           

%                 magData(:,:,nbslices,k) = magDataTmp;
%                 phasData(:,:,nbslices,k) = phasDataTmp;
%                 realData(:,:,nbslices,k) = realDataTmp;
%                 imagData(:,:,nbslices,k) = imagDataTmp;
%                 data(:,:,nbslices,k)=complexDataTmp;  
        
            extra.tVec(k) = d(kk(k)).inversionTime;
            
            s = s+2;
            if (s > nbim)
                break;
            else
            end
                
    end
        
      
    

    
extra.T1Vec = 1:5000; % this range can be reduced if a priori information is available

%data = abs(data);

TI = extra.tVec
save(savename,'data','extra')    

% Display complex data to verfiy the complex conversions:

% disp('display data for the first TI');
% figure;
% imshow(real(data(:,:,1,1)),[]);
% hold;
% 
% disp('display data for the second TI');
% figure;
% imshow(real(data(:,:,1,2)),[]);
% 
% hold;
% 
% disp('display data for the third TI');
% figure;
% imshow(real(data(:,:,1,3)),[]);
% 
% hold;
% 
% disp('display data for the fourth TI');
% figure;
% imshow(real(data(:,:,1,4)),[]);
% 
% hold;

% disp('display d(1).imData');
% figure;
% imshow(d(1).imData(:,:,1,1),[]);



        
        
        
    



