%%
% A MATLAB implementation Psychology and Art Theory-based Features (MH)
% 
%  This set of MATLAB scripts extract the MH features [1] used in [2].
% This code does NOT run itself out-of-the-box since it depends on many
% other programs and the AIC dataset. While we tried our best to make the 
% job as simplified as possible, some of them must be downloaded and 
% compiled separately due to license/compatibility issues. The code to 
% compute Tamura texture features were ported from a C code of [7] under 
% the MIT license included.
%
%  You must cite [1] in order to use these features. In addition, please 
% consider citing [2] if you found this code useful.
%
%  Please let me know if you found any bug(s).
%
% HOW TO USE
%  STEP 1. Download hsy_rgb.zip from http://allan.hanbury.eu/.
%          The file is located under Colour Resources (as of 2015/02/10).
%          Copy the two files into the same directory of this file:
%           hsy2rgb.m, rgb2hsy.m
%
%  STEP 2. Download FastEMD solver from [3]. Compile. Copy following two
%          files into the same directory of this file:
%           emd_hat_gd_metric_mex.mex*, emd_hat_mex.mex* 
%
%  STEP 3. Download Color Descriptor from [4]. Compile. Copy following 
%          files into the same directory of this file:
%           mexColorNaming.mex*, w2c.mat
% 
%  STEP 4. Download color space converter from [5]. Compile. Copy following
%          files into the same directory of this file:
%           colorcalc.mex*, colorspace.m
%
%  STEP 5. Download RGB histogram calculator from [6]. Copy following file
%          into the same directory of this file:
%           rgbhist_fast.m
%
%  STEP 6. You need the MATLAB Computer Vision Toolbox. If not, you need 
%          to obtain a Viola-Johns Face Detector implementation somewhere.
%
%  STEP 7. You need to obtain AIC dataset to properly train the Fuzzy
%          C-Means. While I included pretrained models (AIC_fcm_*.mat), 
%          it is strongly recommended to train them again. Particularly,
%          if you were to reproduce [1], you must re-train it using IAPS
%          dataset instead of the pretrained models. Please check out
%          the codes fuzzyFuncLearn*.m for more details. 
%
% LICENSE
%    A MATLAB implementation Psychology and Art Theory-based Features (MH)
%    Copyright (C) 2014 Sejong Yoon (sjyoon@cs.rutgers.edu)
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% REFERENCE
%  [1] J. Machajdik and A. Hanbury,
%      Affective Image Classification using Features Inspired by
%      Psychology and Art Theory,
%      ACM Intl. Conf. on Multimedia (MM), 2010.
%      URL: http://www.imageemotion.org/
%  [2] S. Yoon and V. Pavlovic, 
%      Sentiment Flow for Video Interestingness Prediction,
%      ACM Intl. Conf. on Multimedia (MM) Workshop (HuEvent), 2014.
%  [3] O. Pele and M. Werman,
%      Fast and Robust Earth Mover's Distances,
%      ICCV 2009.
%      URL: http://www.ariel.ac.il/sites/ofirpele/FastEMD/
%  [4] J. van de Weijer, C. Schmid, J. Verbeek, D. Larlus
%      Learning Color Names for Real-World Applications,
%      IEEE Trans. in Img. Proc. (TIP), vol 18 (7):1512-1524, 2009.
%      URL: http://cat.uab.es/~joost
%  [5] P. Getreuer
%      MATLAB Central File Exchange - Colorspace Transformations
%      URL: .../28790-colorspace-transformations
%  [6] M. K. Reddy,
%      MATLAB Central File Exchange - Color Histogram of an RGB Image.
%      URL: .../43630-color-histogram-of-an-rgb-image
%  [7] T. Minka and R. W. Picard,
%      A Sample Implementation of Tamura Texture Feature in C.
%      URL: http://vismod.media.mit.edu/pub/tpminka/features/
%
% Last Updated: February 10, 2015
%

function fullfeatures = collect_MH_Features(ifname)

%% Initialization
img = imread(ifname);
if ismatrix(img)
    img = repmat(img, [1, 1, 3]);
end
grayimg = rgb2gray(img);
hsyimg = rgb2hsy(img); % code from http://allan.hanbury.eu/

fullfeatures = zeros(114,1);

%% 1-2. color feature (saturation, brightness)
hsyimg_mean = squeeze(mean(mean(hsyimg)));

fullfeatures(1) = hsyimg_mean(2); % Saturation
fullfeatures(2) = hsyimg_mean(3); % Brightness (Luminance)

%% 3-5. pleasure, arousal, dominance
fullfeatures(3) = mean(mean(0.69 * hsyimg(:,:,3) + 0.22 * hsyimg(:,:,2)));
fullfeatures(4) = mean(mean(-0.31 * hsyimg(:,:,3) + 0.60 * hsyimg(:,:,2)));
fullfeatures(5) = mean(mean(0.76 * hsyimg(:,:,3) + 0.32 * hsyimg(:,:,2)));

%% 6-9. hue (Allan Hanbury, Circular statistics applied to colour images, CVWW'03)
H = reshape(hsyimg(:,:,1), [1, size(img,1)*size(img,2)]);
A = sum(cos(H));
B = sum(sin(H));
fullfeatures(6) = atan(A / (B + eps));

R_bar = sqrt(A.^2 + B.^2) / length(H);
fullfeatures(7) = 1 - R_bar;

S = reshape(hsyimg(:,:,2), [1, size(img,1)*size(img,2)]);
A = sum(S .* cos(H));
B = sum(S .* sin(H));
fullfeatures(8) = atan(A / (B + eps));

R_bar = sqrt(A.^2 + B.^2) / (sum(S) + eps);
fullfeatures(9) = 1 - R_bar;

%% 10. colorfulness (Datta, et al. ECCV'06)
fullfeatures(10) = colorfulness_datta(img);

%% 11-21. color names (Joost van de Weijer, et al. CVPR'07)
% black, blue, brown, green, gray, orange, pink, purple, red, white, yellow
fullfeatures(11:21) = colornames_weijer(img);

%% 22-41. Itten color model and 98. level of detail (segmentation)
[fullfeatures(22:41), fullfeatures(98)] = itten(img);

%% 42-70. Wang features with area statistics
fullfeatures(42:70) = wang(img);

%% 71-73. Tamura texture (coarseness, contrast, directionality)
tamFeatures = tamuraFeatures(grayimg);
fullfeatures(71:73) = tamFeatures(1:3);

%% 74-85. wavelet texture and 99-101. low depth of field
[fullfeatures(74:85), fullfeatures(99:101)] = waveletFeatures(hsyimg);

%% 86-97. GLCM (contrast, correlation, energy, homogeneity)
glcmH = graycoprops(graycomatrix(hsyimg(:,:,1)));
glcmS = graycoprops(graycomatrix(hsyimg(:,:,2)));
glcmY = graycoprops(graycomatrix(hsyimg(:,:,3)));
fullfeatures(86) = glcmH.Contrast;
fullfeatures(87) = glcmS.Contrast;
fullfeatures(88) = glcmY.Contrast;
fullfeatures(89) = glcmH.Correlation;
fullfeatures(90) = glcmS.Correlation;
fullfeatures(91) = glcmY.Correlation;
fullfeatures(92) = glcmH.Energy;
fullfeatures(93) = glcmS.Energy;
fullfeatures(94) = glcmY.Energy;
fullfeatures(95) = glcmH.Homogeneity;
fullfeatures(96) = glcmS.Homogeneity;
fullfeatures(97) = glcmY.Homogeneity;

%% 102-107. dynamics
BW = edge(grayimg, 'canny');
[H, theta, rho] = hough(BW);
P = houghpeaks(H, 5, 'threshold', ceil(0.3*max(H(:))));
lines = houghlines(BW, theta, rho, P, 'FillGap', 5, 'MinLength', 7);
% subplot(121); imshow(BW);
% subplot(122); 
% hold on;
% max_len = 0;
static_abs_deg_sum = 0;
slant_abs_deg_sum = 0;
static_rel_deg_sum = 0;
slant_rel_deg_sum = 0;
static_len_sum = 0;
slant_len_sum = 0;
static_count = 0;
slant_count = 0;
for k = 1 : length(lines)
%     xy = [lines(k).point1; lines(k).point2];
%     plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% 
%     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    len = norm(lines(k).point1 - lines(k).point2);
    t = lines(k).theta;
    if (t > -15 && t < 15)  || (t > 75 || t < -75)
        static_len_sum = static_len_sum + len;
        static_abs_deg_sum = static_abs_deg_sum + t;
        static_rel_deg_sum = static_rel_deg_sum + (t + 90) / 180;
        static_count = static_count + 1;
    else
        slant_len_sum = slant_len_sum + len;
        slant_abs_deg_sum = slant_abs_deg_sum + t;
        slant_rel_deg_sum = slant_rel_deg_sum + (t + 90) / 180;
        slant_count = slant_count + 1;
    end
%     if len > max_len
%         max_len = len;
%         xy_long = xy;
%     end
%     fprintf('%d\n', t);
%     pause;
end
% plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','blue');
% hold off;
fullfeatures(102:104) = [static_abs_deg_sum, static_rel_deg_sum, static_len_sum];
fullfeatures(105:107) = [slant_abs_deg_sum, slant_rel_deg_sum, slant_len_sum];
fullfeatures(102:104) = fullfeatures(102:104) ./ (static_count+eps);
fullfeatures(105:107) = fullfeatures(105:107) ./ (slant_count+eps);

%% 108-110. rule of thirds
[h, w, c] = size(hsyimg);
hlb = floor(h/3);
hub = floor(2*h/3);
wlb = floor(w/3);
wub = floor(2*w/3);
XY = ((hub-hlb)*(wub-wlb));
fullfeatures(108) = sum(sum(hsyimg(hlb:hub,wlb:wub,1))) * 9 / XY;
fullfeatures(109) = sum(sum(hsyimg(hlb:hub,wlb:wub,2))) * 9 / XY;
fullfeatures(110) = sum(sum(hsyimg(hlb:hub,wlb:wub,3))) * 9 / XY;

%% 111-112. faces (Require: Computer Vision Toolbox)
faceDetector = vision.CascadeObjectDetector();
bbox = step(faceDetector, img);

if size(bbox,1) > 0
    fullfeatures(111) = size(bbox,1);
    fullfeatures(112) = max(bbox(:,3) .* bbox(:,4)) / (h*w);
end

%% 113-114. skin
fullfeatures(113:114) = skinFeatures(img, bbox);

end
