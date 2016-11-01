function [features, nSeg] = itten(img)
% ITTEN
%
%   Computes Itten color model features
%
% INPUT
%   img      : image to extract features
%
% REFERENCES
%
%   J. Itten and E. Haagen, The Art of Color: The Subjective Experience and
%   Objective Rationale of Color, 1974.

features = zeros(20, 1);

% Do watershed segmentation and calculate the needed statistics
[segLab, segAvg, segVar] = watershedSegment(img);

nSeg = size(segAvg,1);
segs = unique(segLab);
segSize = zeros(length(segs),1);
for ids = 1 : length(segs)
    segSize(ids) = sum(segLab == segs(ids));
end
segSize = segSize ./ sum(segSize);

% load membership function
fcm_model = fuzzyFuncLearnItten();

% predict membership function
[Us, Uy] = fuzzyFuncPredict(fcm_model, segAvg);

% Calculate features

% (1) Contrast of Brightness
yContrast = sqrt(var(segSize' * Uy,1,2));
features(1) = mean(yContrast);

% (2) Contrast of Saturation
sContrast = sqrt(var(segSize' * Us,1,2));
features(2) = mean(sContrast);

% (3) Contrast of Hue
hueDist = pdist(segAvg(:,1));
hueDist = min(hueDist, 2*pi - hueDist);
hContrast = 2 * sin(hueDist / 2);
features(3) = mean(hContrast);

% (4) Contrast of complements
features(4) = mean(hueDist);

% (5) Contrast of Warmth
warms = (0 <= segAvg(:,1) & segAvg(:,1) < deg2rad(140)) | ...
    (deg2rad(320) <= segAvg(:,1) & segAvg(:,1) < 2*pi);
colds = (deg2rad(140) <= segAvg(:,1) & segAvg(:,1) < deg2rad(320));
MuWarmth = zeros(size(segAvg,1),3);
MuWarmth(colds,1) = cos(segAvg(colds,1) - deg2rad(230));
MuWarmth(warms,3) = cos(segAvg(warms,1) - deg2rad(50));
MuWarmth(:,2) = 1 - (MuWarmth(:,3) + MuWarmth(:,1));

numer = MuWarmth * MuWarmth';
denom = sum(MuWarmth .^ 2, 2);
denom = sqrt(denom * denom');
WarmthContrast = numer ./ denom;

features(5) = mean(mean(WarmthContrast));

% (6) Harmony
binEdges = deg2rad([0:30:360] + 15); % bracket is *necessary*!
hueCount = histc(segAvg(:,1) + deg2rad(15), binEdges);
nHueCount = hueCount ./ sum(hueCount);
nHueCount(nHueCount < 0.05) = 0; % drop bins with less than 5% support

n = length(find(nHueCount > 0));
regularIntAng = pi * (n - 2) / n;

% calculate angle differences
hueAngles = binEdges(2:end);
diffHue = zeros(n,1);
switch n
    case {1,2}
        diffHue = regularIntAng;
    otherwise,
        % convert polar to Cartesian
        xy = [cos(hueAngles'), sin(hueAngles')];
        % pick only found
        xy = xy(nHueCount > 0,:);
        % append first and last to the last and first
        xy = [xy(end,:); xy; xy(1,:)];
        % get angles between the two lines
        for idp = 2 : (n+1)
            org = xy(idp,:);
            p1 = xy(idp-1,:) - org;
            p2 = xy(idp+1,:) - org;
            dxy = p2 - p1;
            diffHue(idp-1) = abs(regularIntAng - (atan2(dxy(2), dxy(1)) + pi));
        end
end
features(6) = mean(diffHue);

% (7) Hue count
features(7) = mean(hueCount);

% (8) Hue spread
features(8) = mean(segVar(:,1));

% (9) Total area of warm
features(9) = mean(segSize(warms));

% (10) Total area of cold
features(10) = mean(segSize(colds));

% Max of each

% (11) Contrast of Brightness
features(11) = max(yContrast); % since yContrast is a scalar, it's the same 

% (12) Contrast of Saturation
features(12) = max(sContrast); % since sContrast is a scalar, it's the same

% avoid some dangerous cases (blank frame)
if isnan(features(3))
    features(isnan(features)) = 0;
    return;
end

% (13) Contrast of Hue
features(13) = max(hContrast);

% (14) Contrast of complements
features(14) = max(hueDist);

% (15) Contrast of Warmth (NOTE: This will be almost always 1...)
features(15) = max(max(WarmthContrast));

% (16) Harmony
features(16) = max(diffHue);

% (17) Hue count
features(17) = max(hueCount);

% (18) Hue spread
features(18) = max(segVar(:,1));

% (19) Area of warm
if sum(warms) > 0
    features(19) = max(segSize(warms));
end

% (20) Area of cold
if sum(warms) > 0
    features(20) = max(segSize(colds));
end

end
