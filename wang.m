function features = wang(img)
% WANG
%
%   Computes fuzzy histogram features proposed by Wang, et al.
%
% INPUT
%   img      : image to extract features
%
% REFERENCES
%   W. Wang, et al., Image Retrieval by Emotional Semantics: A Study of
%   Emotional Space and Feature Extraction, IEEE Intl. Conf. on Sys. Man,
%   and Cybernetics, 2006.

features = zeros(29,1);

fcm_model = fuzzyFuncLearnWang();

if ismatrix(img)
    img = repamt(img, [1, 1, 3]);
end

[h, w, c] = size(img);

totalPixels = h*w;

grayimg = double(rgb2gray(img))./255;

colorTransform = makecform('srgb2lab');
labimg = applycform(img, colorTransform);
labimg = reshape(labimg, [h*w, c]);

img = reshape(rgb2hsy(img), [h*w, c]);

[~, Uy] = fuzzyFuncPredict(fcm_model, img);

% Factor 1 (10)
[~, Uy] = max(Uy,[],2);
veryDark = (Uy == 1);
dark = (Uy == 2);
middle = (Uy == 3);
light = (Uy == 4);
veryLight = (Uy == 5);

warms = (0 <= img(:,1) & img(:,1) < deg2rad(140)) | ...
    (deg2rad(320) <= img(:,1) & img(:,1) < 2*pi);
cool = (deg2rad(140) <= img(:,1) & img(:,1) < deg2rad(320));

features(1) = sum(veryDark & warms);
features(2) = sum(veryDark & cool);
features(3) = sum(dark & warms);
features(4) = sum(dark & cool);
features(5) = sum(middle & warms);
features(6) = sum(middle & cool);
features(7) = sum(light & warms);
features(8) = sum(light & cool);
features(9) = sum(veryLight & warms);
features(10) = sum(veryLight & cool);

% Factor 2 (7)
areaS1 = img(:,2) < (10/100);
areaS2 = img(:,2) >= (20/100) & img(:,2) <= (27/100);
areaS3 = img(:,2) >= (10/100) & img(:,2) < (27/100);
areaS4 = img(:,2) >= (27/100) & img(:,2) < (51/100);
areaS5 = img(:,2) >= (27/100) & img(:,2) <= (51/100);
areaS6 = img(:,2) > (51/100);

LowS = areaS1 + areaS2 .* (0.27 - img(:,2)) ./ 17;
MiddleS = areaS3 .* (img(:,2) - 0.1) ./ 17 + areaS4 .* (0.51 - img(:,2)) ./ 24;
HighS = areaS5 .* (img(:,2) - 0.27) ./ 24 + areaS6;

features(11) = sum(warms & HighS);
features(12) = sum(cool & HighS);
features(13) = sum(warms & MiddleS);
features(14) = sum(cool & MiddleS);
features(15) = sum(warms & LowS);
features(16) = sum(cool & LowS);

meanLabImg = mean(labimg(:,2:3), 1);
meanLabImg = repmat(meanLabImg, [h*w, 1]);
lab_mean = sum(sum((double(labimg(:,2:3)) - meanLabImg).^2));
features(17) = sqrt(lab_mean / (h*w-1));

% Factor 3 (2)
features(18) = var(img(:,3)); % brightness contrast

sf = fspecial('sobel');
gih = double(imfilter(grayimg, sf));
giv = double(imfilter(grayimg, sf'));
edgeimg = sort(reshape(sqrt(giv.^2+gih.^2), [h*w 1]), 1, 'descend');
edgeimg = edgeimg(1:ceil(h*w*0.005)); % sample only top 0.5%
features(19) = mean(edgeimg); % sharpness

% Noramlize!
features(1:18) = features(1:18) ./ totalPixels;

% Area statistics
features(20) = sum(veryDark) / totalPixels;
features(21) = sum(dark) / totalPixels;
features(22) = sum(middle) / totalPixels;
features(23) = sum(light) / totalPixels;
features(24) = sum(veryLight) / totalPixels;
features(25) = sum(warms) / totalPixels;
features(26) = sum(cool) / totalPixels;
features(27) = sum(HighS) / totalPixels;
features(28) = sum(MiddleS) / totalPixels;
features(29) = sum(LowS) / totalPixels;

end
