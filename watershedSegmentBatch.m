function [segLabel, segAvg, segVar] = watershedSegmentBatch(imgpath, ext)
%   WATERSHED_SEGMENT_BATCH
%
%   Watershed segmentation on all images with extension ext in impath.

files = dir([imgpath '/*.' ext]);

nFiles = length(files);
segLabel = cell(nFiles,1);
segAvg = cell(nFiles,1);
segVar = cell(nFiles,1);

for idf = 1 : nFiles
    img = imread([imgpath '/' files(idf).name]);
    if ismatrix(img)
        img = repmat(img, [1 1 3]);
    end
    
    L = watershedSegmentBatchOne(img);
    [h, w] = size(L);
    L = reshape(L, [h*w 1]);
    
    segLabel{idf} = L;

    hsyimg = rgb2hsy(img);
    hsyimg = reshape(hsyimg, [h*w 3]);

    segs = unique(L);
    segAvg{idf} = zeros(length(segs),3);
    segVar{idf} = zeros(length(segs),3);
    for ids = 1 : length(segs)
        meanhsy = mean(hsyimg(L == segs(ids),:), 1);
        varhsy = var(hsyimg(L == segs(ids),:), [], 1);
        segAvg{idf}(ids,:) = meanhsy; % mean HSL for each cluster
        segVar{idf}(ids,:) = varhsy;
    end
    
    fprintf('Processed %d / %d files.\n', idf, nFiles);
end

end

function L = watershedSegmentBatchOne(img)
% Follows the MATLAB example of watershed segmentation

% watershed segmentation
gimg = rgb2gray(img);

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(gimg), hy, 'replicate');
Ix = imfilter(double(gimg), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);

se = strel('disk', 3);
Iobr = imreconstruct(imerode(gimg, se), gimg);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);

se2 = strel(ones(3,3));
fgm = bwareaopen(imerode(imclose(imregionalmax(Iobrcbr), se2), se2), 3);

DL = watershed(bwdist(im2bw(Iobrcbr, graythresh(Iobrcbr))));
bgm = DL == 0;

gradmag2 = imimposemin(gradmag, bgm | fgm);
L = watershed(gradmag2);
% Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
% imshow(Lrgb);

end