function [features] = skinFeatures(img, bbox)
%   SKIN_FEATURES
%
%   If bbox is not empty (i.e. face detected by Viola-Johns), use
%   simplified version of [1] by [2]. Otherwise, use [3]
%
% REFERENCES
%   [1] C. Liensberger, et al., Color-Based and Context-Aware Skin 
%   Detection for Online Video Annotation, MMSP'09.
%   [2] J. Machajdik and A. Hanbury, Affective Image Classification using
%   Features Inspired by Psychology and Art Theory, ACM MM'10.
%   [3] D. Chai and K.N. Ngan, Locating Facial Region of a 
%   Head-and-Shoulders Color Image, In Proc. Of IEEE Region Ten Conference, 
%   vol. 2, 421- 4124, 1999.

features = zeros(2,1);

[h, w, ~] = size(img);

ycbcrimg = rgb2ycbcr(img);
cb = ycbcrimg(:,:,2);
cr = ycbcrimg(:,:,3);

if size(bbox,1) == 0
    cond1 = cb >= 77 & cb <= 127;
    cond2 = cr >= 133 & cr <= 173;
    
    skinPixels = cond1 & cond2;        
else
    maxCb = [];
    minCb = [];
    maxCr = [];
    minCr = [];
    for idf = 1 : size(bbox,1)
        hrange = bbox(idf,2):(bbox(idf,2)+bbox(idf,4));
        wrange = bbox(idf,1):(bbox(idf,1)+bbox(idf,3));
        nh = length(hrange);
        nw = length(wrange);
        maxCb = [maxCb; reshape(ycbcrimg(hrange, wrange, 2), [nh*nw, 1])];
        minCb = [minCb; reshape(ycbcrimg(hrange, wrange, 2), [nh*nw, 1])];
        maxCr = [maxCr; reshape(ycbcrimg(hrange, wrange, 3), [nh*nw, 1])];
        minCr = [minCr; reshape(ycbcrimg(hrange, wrange, 3), [nh*nw, 1])];
    end
    nCb = round(length(maxCb) * 0.3);
    nCr = round(length(maxCb) * 0.175);
    maxCb = sort(maxCb,1,'descend');
    minCb = sort(minCb,1,'ascend');
    maxCr = sort(maxCr,1,'descend');
    minCr = sort(minCr,1,'ascend');
    maxCb = mean(maxCb(1:nCb));
    minCb = mean(minCb(1:nCb));
    maxCr = mean(maxCr(1:nCr));
    minCr = mean(minCr(1:nCr));
    
    skinPixels = (cb >= minCb & cb <= maxCb & cr >= minCr & cr <= maxCr);
end

features(1) = sum(sum(skinPixels)) / (h*w); % normalized size of skin area

if size(bbox,1) > 0
    facearea = bbox(:,3)' * bbox(:,4);
    features(2) = features(1) ./ facearea; % ratio of non-facial skin to face
else
    features(2) = 1.0; % assume all pixels are faces
end

end
