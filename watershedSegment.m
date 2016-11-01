function [segLabel, segAvg, segVar] = watershedSegment(img)
%   WATERSHED_SEGMENT
%
%   Watershed segmentation on given image...just following MATLAB example

if ismatrix(img)
    img = repmat(img, [1 1 3]);
end

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

[h, w] = size(L);
segLabel = reshape(L, [h*w 1]);

hsyimg = rgb2hsy(img);
hsyimg = reshape(hsyimg, [h*w 3]);

segs = unique(segLabel);
segAvg = zeros(length(segs),3);
segVar = zeros(length(segs),3);
for ids = 1 : length(segs)
    meanhsy = mean(hsyimg(segLabel == segs(ids),:), 1);
    varhsy = var(hsyimg(L == segs(ids),:), [], 1);
    segAvg(ids,:) = meanhsy; % mean HSL for each cluster
    segVar(ids,:) = varhsy;
end
    
end
