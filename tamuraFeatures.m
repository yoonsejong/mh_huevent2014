function features = tamuraFeatures(img, param)
%   TAMURA_FEATURES
%
%   Calculates the three Tamura texture featuress [1] for a given grayscale 
%   image. This is a MATLAB port of C code from [2].
%
% INPUT
%   img          : grayscale image matrix
%   param.       : (optional) structured parameter
%    max_scale   : maximum level for coarseness (def:5)
%    weight      : weight for coarseness (def:1)
%    histo_bins  : angle histogram bins for directionality (def:16)
%    histo_thres : edge magnitude threshold for directionality (def:12)
%    peak_thres  : peak threshold for directionality (def:2)
%
% OUPUT
%   features     : four dimensional vector each for coarseness, contrast,
%                  directionality and orientation of Tamura features
%
%  REFERENCE
%   [1] H. Tamura, Textural Features Corresponding to Visual Perception,
%   IEEE Trans. on Sys. Man, AND Cybernetics, Vol. SMC-8, No. 6, June 1978.
%   [2] http://vismod.media.mit.edu/pub/tpminka/features/

if nargin < 1
    % Output should be [0.3292, 0.4394, 0.2631, 1.0000]
    load('tamuratest.mat');
    img = tamuratest;
end
if nargin < 2
    max_scale = 5;
    weight = 1;
    histo_bins = 16;
    histo_thres = 12;
    peak_thres = 2;
else
    if isfield(param, 'max_scale'), max_scale = param.max_scale; end;
    if isfield(param, 'weight'), weight = param.weight; end;
    if isfield(param, 'histo_bins'), histo_bins = param.histo_bins; end;
    if isfield(param, 'histo_thres'), histo_thres = param.histo_thres; end;
    if isfield(param, 'peak_thres'), peak_thres = param.peak_thres; end;
end

if ~ismatrix(img)
    error('Sorry, only grayscale images!');
end

features = zeros(4, 1);

% calculate
features(1) = tamuraCoarseness(img, max_scale, weight);
features(2) = tamuraContrast(img);
features(3:4) = tamuraDirectionality(img, histo_bins, histo_thres, peak_thres);

end

%% Coarseness
function feature = tamuraCoarseness(img, max_scale, weight)
for scale = 1 : max_scale
    img = nbr_avg(img, scale);
    diff = nbr_diff(img, scale);
    diff = diff .* (weight.^(scale-1));
    if ~exist('best', 'var')
        best = diff;
        scales = ones(size(diff)) * 2.^(scale);
    else
        step = floor(2.^(scale-1));
        [diffH, diffW] = size(diff);
        best = best(step:(diffH+step),step:(diffW+step));
        scales = scales(step:(diffH+step),step:(diffW+step));
        
        % NOTE: the sizes of best and diff are different!
        for y = 1 : diffH
            for x = 1 : diffW
                if best(y,x) <= diff(y,x)
                    best(y,x) = diff(y,x);
                    scales(y,x) = 2.^(scale);
                end
            end
        end
    end
end
crs = mean(mean(scales));
feature = crs / 2.^(max_scale);
end

function res = nbr_diff(img, factor)
step = floor(2.^factor);
% [height, width] = size(img);
% res1 = zeros(height-step, width-step);
% for y = 1 : (height-step)
%     for x = 1 : (width-step)
%         h_diff = img(y,x) - img(y,x+step);
%         v_diff = img(y,x) - img(y+step,x);
%         res1(y,x) = max(abs(h_diff), abs(v_diff));
%     end
% end

% Long live MATLAB!
FH = zeros(step*2+1,step*2+1);
FH(step+1,step+1) = 1;
FH(step+1,1) = -1;
FV = zeros(step*2+1,step*2+1);
FV(step+1,step+1) = 1;
FV(1,step+1) = -1;
fimh = imfilter(img, FH, 'conv'); 
fimv = imfilter(img, FV, 'conv');
res = max(abs(fimh), abs(fimv));
res = res(1:(end-step),1:(end-step));
end

function res = nbr_avg(img, factor)
step = floor(2.^(factor-1));
% [height, width] = size(img);
% res1 = zeros(height-step,width-step);
% for y = 1 : (height-step)
%     for x = 1 : (width-step)
%         res1(y,x) = (img(y,x) + img(y,x+step) + ...
%             img(y+step,x) + img(y+step,x+step)) / 4;
%     end
% end

% Long live MATLAB!
F = zeros(step*2+1,step*2+1);
F(1,1) = 1;
F(step+1,1) = 1;
F(1,step+1) = 1;
F(step+1,step+1) = 1;
res = double(imfilter(img, F, 'conv') ./ 4); 
res = res(1:(end-step),1:(end-step));
end

%% Contrast
function feature = tamuraContrast(img)
img_mean = mean(mean(img));
[height, width] = size(img);
img_diff = double(img) - img_mean .* ones(size(img));
img_var = img_diff .* img_diff;
img_norm4 = img_var .* img_var;
img_var = sum(sum(img_var)) / (height*width);
img_norm4 = sum(sum(img_norm4)) / (height*width);
if img_norm4 == 0
    feature = 0;
else
    % normalize by maximum value
    feature = img_var / img_norm4.^(0.25) / 128;
end
end

%% Directionality
function features = tamuraDirectionality(img, histo_bins, histo_thres, peak_thres)
features = zeros(2,1);

% compute direction histogram
histogram = dir_histogram(img, histo_bins, histo_thres);

% find extrema of histogram
extrema = zeros(histo_bins,1);
is_peak = zeros(histo_bins,1);
n_extrema = 0;

% compute difference vector
diff = histogram(2:end) - histogram(1:(end-1));
diff = [diff; (histogram(1) - histogram(end))];
if (diff(histo_bins) * diff(1) < 0) || ...
        (diff(histo_bins) * diff(1) == 0 && diff(histo_bins) + diff(1) ~= 0)
    n_extrema = n_extrema + 1;
    extrema(n_extrema) = 1;
    is_peak(n_extrema) = diff(1) < 0;
end
for i = 1 : (histo_bins-1)
    if (diff(i) * diff(i+1) < 0) || ...
            (diff(i) * diff(i+1) == 0 && diff(i) + diff(i+1) ~= 0)
        n_extrema = n_extrema + 1;
        extrema(n_extrema) = i+1;
        is_peak(n_extrema) = diff(i+1) < 0;
    end
end

% extract salient peaks
pv = zeros(n_extrema,3);
j = 0;
num_peaks = 0;
for i = 1 : n_extrema
    if is_peak(i)
        if j == 0
            k = n_extrema;
            while is_peak(k)
                k = k -1;
            end
            pv(num_peaks+1,1) = extrema(k);
            j = j + 1;
        end
        if j == 1
            pv(num_peaks+1,2) = extrema(i);
            j = j + 1;
        end
    else
        if j < 2
            pv(num_peaks+1,1) = extrema(i);
            j = 1;
        else
            pv(num_peaks+1,3) = extrema(i);
            if valid_peak(pv(num_peaks+1,:), histogram, peak_thres)
                num_peaks = num_peaks + 1;
            end
            pv(num_peaks+1,1) = extrema(i);
            j = 1;
        end
    end
end
if j == 2
    k = 1;
    while is_peak(k)
        k = k + 1;
    end
    pv(num_peaks+1,3) = extrema(k);
end
if j < 2 || ~valid_peak(pv(num_peaks+1,:), histogram, peak_thres)
    num_peaks = num_peaks - 1;
end
num_peaks = num_peaks + 1;

features(2) = num_peaks;

% compute 2nd moment about peaks
if num_peaks == 0
    features(1) = 1;
    return;
end
dir = 0;
for i = 0 : (num_peaks-1)
    dir = dir + sd_histo(histogram, histo_bins, pv(i+1,1), pv(i+1,3));
end

features(1) = dir / 5; % normalize by maximum value
end

function histo = dir_histogram(img, histo_bins, threshold)
% TODO: We should be able to utilize edge() and histc() here.
histo = zeros(histo_bins,1);
sumpix = 0;
dont_care = 0;
[height, width] = size(img);
for y = 2 : (height-1)
    for x = 2 : (width-1)
        delh = (img(y-1,x+1) + img(y,x+1) + img(y+1,x+1)) - ...
            (img(y-1,x-1) + img(y,x-1) + img(y+1,x-1));
        delv = (img(y-1,x-1) + img(y-1,x) + img(y-1,x+1)) - ...
            (img(y+1,x-1) + img(y+1,x) + img(y+1,x+1));
        delG = (abs(delh) + abs(delv)) / 2;    
        if delG >= threshold
            theta = atan2(double(delv), double(delh));
            if theta < 0, 
                theta = theta + pi;
            elseif theta >= pi, 
                theta = theta - pi;
            end
            bin = floor(theta * histo_bins / pi + 0.5) + 1;
            if bin == (histo_bins+1), bin = 1; end;
            if bin < 1 || bin > histo_bins
                error('Bin error %d', bin);
            end
            histo(bin) = histo(bin) + 1;
            sumpix = sumpix + 1;
        else
            dont_care = dont_care + 1;
        end
    end
end

sumpix = sumpix + dont_care;
dont_care = dont_care / histo_bins;
if sumpix
    histo = (histo + dont_care) ./ sumpix;
end
end

function res = valid_peak(peak, histogram, peak_thres)
max_valley = histogram(peak(1));
if histogram(peak(3)) > max_valley
    max_valley = histogram(peak(3));
end
sharpness = histogram(peak(2)) / max_valley;
res = (sharpness > peak_thres);
end

function res = sd_histo(histogram, histo_bins, start, stop)
length = stop - start;
if stop <= start
    length = length + histo_bins; 
end
interval = zeros(length,1);

% copy the sub-histogram into the interval
if stop <= start
    for i = (start+1) : histo_bins
        interval(i - start) = histogram(i);
    end
    for i = 1 : stop
        interval(histo_bins - start + i) = histogram(i);
    end
else
    for i = (start+1) : stop
        interval(i - start) = histogram(i);
    end
end

% normalize
interval = interval ./ sum(interval);

% find mean
mp = ceil(sum(interval .* [1:length]') - 0.5);

% compute variance
res = sum(interval .* ([1:length]' - mp).^2);

res = sqrt(res);
end
