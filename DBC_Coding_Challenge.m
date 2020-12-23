%% Dynamic Brain Circuits coding challenge




%% First clone the seqNMF repository from github: 
% https://github.com/elifesciences-publications/seqNMF.git


%% Load the data

clc
cd F:\CL\20190807173901\S1

I = h5read('image_stream.hdf5','/raw_images');





%% Draw the ROI
% Large changes in pixel values can often be observed over the central
% sinus or at the edges of the transcranial window. Draw an ROI that masks
% both hemispheres of cortex but excludes these problematic areas. 

X = squeeze(I(2,:,:,:));
X = imrotate(X(:,:,1000:1499),270);
X = imresize(X,0.25);
roi = draw_roi(X,2);
X = double(X).*roi;


%% Calculate the fractional change in fluorescence (dF/F0) for all pixels
% This can be calculated as (F-F0)/F0. Where F is the vector of pixel 
% values, and F0 is the "baseline" pixel value (mean across time). 
%
% A better approach uses a moving baseline pixel value (10s). This can help
% account for artifacts such as photobleaching.


% 9.75s baseline for dFF
fs = 15; % frame rate
baseline = movmean(X,round(9.75*fs),3);
dff = double((X-baseline)./baseline);    
dff = zscore(dff,[],3);


%% Filter the data with a bandpass filter (0.1-4Hz). 

[b,a] = butter(1, [0.1 4]/(fs/2));
fX = dff;
fX(isnan(fX)) = 0;
fX = filtfilt(b,a,permute(fX,[3 1 2]));
fX2 = permute(fX,[2 3 1]);

%% Format data to be NxT
% Data should be 4096x500
% ROI should be a binary matrix of size 4096x500

X2 = reshape(fX2,size(fX2,1)*size(fX2,2),[]);    
M = reshape(roi, size(roi,1)*size(roi,2),[]);

%% Normalize data between 0 and 1
% Traditional NMF approaches require that the matrix to be factorized has 
% no negative elements.

X2 = (X2 - min(X2(:))) ./ (max(X2(:))-min(X2(:)));


%% Run seqNMF
% Input arguments have been adjusted for the data.
% Do not change

tic
[W, H, cost, loadings, power] = seqNMF(X2, ...
    'K',15, 'L',15, 'lambda',0.00005, 'showPlot',1, 'maxiter',100, ...
    'tolerance',0, 'shift',0, 'lambdaL1W',0, 'lambdaL1H',1, ...
    'W_fixed',0, 'W_init',nan, 'H_init',nan, 'SortFactors',0, ...
    'lambdaOrthoH',1, 'lambdaOrthoW',0, 'useWupdate',1, 'M',M);
toc

%% Visualize outputs


motif = cell(1, size(W,2));

for i = 1:size(W,2)  

    tmp = squeeze(W(:,i,:));
    motif{i} = reshape(tmp, 64, 64, []);
    motif{i} = (motif{i}-min(motif{i}(:))) ./ (max(motif{i}(:)) - min(motif{i}(:)));

end

row1 = [];
row2 = [];
row3 = [];

sz = length(motif)/3;
for i = 1:sz

    row1 = cat(2,row1, motif{i});
    row2 = cat(2,row2, motif{sz+i});
    row3 = cat(2,row3, motif{2*sz+i});
end

tiledMotifs = imgaussfilt3(cat(1,row1,row2,row3), [1, 1, 0.1]);
tiledMotifs(isnan(tiledMotifs)) = 0;
cmap = jet(256);
figure, sliceViewer(tiledMotifs, 'Colormap', cmap)


%%
function roi = draw_roi(img,numroi)
% Draw regions of interest on image data. This script was adapted from 
% https://www.mathworks.com/matlabcentral/answers/uploaded_files/113201/draw_multiple_polygons.m
% 
% Inputs:
%   img          (required, image data - dimensions: Height x Width x Time)
%   numroi       (optional, number of ROIs - default = 1)
%
% Outputs:
%   roi          (roi mask)
% 
% Usage:  roi = draw_roi(img, numroi);

if nargin < 2 || isempty(numroi), numroi = 2; end

roihandle = figure();
subplot(1,2,1), imagesc(mean(img,3)), title('Original')
subplot(1,2,2), imagesc(mean(img,3)), title('ROIs overlaid'), 
colormap gray, hold on

% Ask user to draw freehand mask.
message = sprintf(['Left click to draw vertices in the left image.',...
    '\nRight click to finish. Vertices can added by pressing A.',...
    '\nDouble click in the middle to accept.']);

answer = questdlg(message,'','Quit','Continue','Continue');

switch answer
    case 'Continue'
        regionCount = 0;
        roi = false(size(img));
        while regionCount < numroi
            regionCount = regionCount + 1;

            subplot(1, 2, 1); % Switch to img axes.
            [single_roi, xi, yi] = roipoly();

            % Draw the polygon over the img on the right.
            subplot(1, 2, 2);
            plot(xi, yi, 'r-', 'LineWidth', 2);

            % Create the combined roi mask
            roi = roi | single_roi;
        end
        roi = single(roi);
        close(roihandle)
    case 'Quit'
        roi = ones(size(img),'single');
end
end
