function [wfeatures, ldoffeatures] = waveletFeatures(hsyimg)
%   WAVELET_FEATURES
%
%   Computes three-level Daubechies wavelet features proposed by Datta, et
%   al. in ECCV'10

%% Wavelet features
wfeatures = zeros(12,1);

imgH = hsyimg(:,:,1);
imgS = hsyimg(:,:,2);
imgY = hsyimg(:,:,3);

% level 1
[LL1H, HL1H, LH1H, HH1H] = dwt2(imgH, 'db1');
[LL1S, HL1S, LH1S, HH1S] = dwt2(imgS, 'db1');
[LL1Y, HL1Y, LH1Y, HH1Y] = dwt2(imgY, 'db1');

numH = sum(sum(HL1H)) + sum(sum(LH1H)) + sum(sum(HH1H));
numS = sum(sum(HL1S)) + sum(sum(LH1S)) + sum(sum(HH1S));
numY = sum(sum(HL1Y)) + sum(sum(LH1Y)) + sum(sum(HH1Y));

[r, c] = size(HL1H);
S = r*c*3;

wfeatures(1) = numH / S;
wfeatures(2) = numS / S;
wfeatures(3) = numY / S;

imgH = [LL1H, HL1H; LH1H, HH1H];
imgS = [LL1S, HL1S; LH1S, HH1S];
imgY = [LL1Y, HL1Y; LH1Y, HH1Y];

% level 2
[LL2H, HL2H, LH2H, HH2H] = dwt2(imgH, 'db1');
[LL2S, HL2S, LH2S, HH2S] = dwt2(imgS, 'db1');
[LL2Y, HL2Y, LH2Y, HH2Y] = dwt2(imgY, 'db1');

numH = sum(sum(HL2H)) + sum(sum(LH2H)) + sum(sum(HH2H));
numS = sum(sum(HL2S)) + sum(sum(LH2S)) + sum(sum(HH2S));
numY = sum(sum(HL2Y)) + sum(sum(LH2Y)) + sum(sum(HH2Y));

[r, c] = size(HL2H);
S = r*c*3;

wfeatures(4) = numH / S;
wfeatures(5) = numS / S;
wfeatures(6) = numY / S;

imgH = [LL2H, HL2H; LH2H, HH2H];
imgS = [LL2S, HL2S; LH2S, HH2S];
imgY = [LL2Y, HL2Y; LH2Y, HH2Y];

% level 3
[~, HL3H, LH3H, HH3H] = dwt2(imgH, 'db1');
[~, HL3S, LH3S, HH3S] = dwt2(imgS, 'db1');
[~, HL3Y, LH3Y, HH3Y] = dwt2(imgY, 'db1');

numH = sum(sum(HL3H)) + sum(sum(LH3H)) + sum(sum(HH3H));
numS = sum(sum(HL3S)) + sum(sum(LH3S)) + sum(sum(HH3S));
numY = sum(sum(HL3Y)) + sum(sum(LH3Y)) + sum(sum(HH3Y));

[r, c] = size(HL3H);
S = r*c*3;

wfeatures(7) = numH / S;
wfeatures(8) = numS / S;
wfeatures(9) = numY / S;

wfeatures(10) = sum(wfeatures(1:3));
wfeatures(11) = sum(wfeatures(4:6));
wfeatures(12) = sum(wfeatures(7:9));

%% Low Depth of Field
ldoffeatures = zeros(3,1);

rlb = floor(r / 4);
rub = floor(r * 3 / 4);

clb = floor(c / 4);
cub = floor(c * 3 / 4);

numHi = sum(sum(HL3H(rlb:rub,clb:cub))) ...
    + sum(sum(LH3H(rlb:rub,clb:cub))) ...
    + sum(sum(HH3H(rlb:rub,clb:cub)));
numSi = sum(sum(HL3S(rlb:rub,clb:cub))) ...
    + sum(sum(LH3S(rlb:rub,clb:cub))) ...
    + sum(sum(HH3S(rlb:rub,clb:cub)));
numYi = sum(sum(HL3Y(rlb:rub,clb:cub))) ...
    + sum(sum(LH3Y(rlb:rub,clb:cub))) ...
    + sum(sum(HH3Y(rlb:rub,clb:cub)));

ldoffeatures(1) = numHi / numH;
ldoffeatures(2) = numSi / numS;
ldoffeatures(3) = numYi / numY;

end
