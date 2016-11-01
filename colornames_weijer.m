function out = colornames_weijer(img)
% COLORNAMES_WEIJER
%
%   Wrapper function to use color names proposed by Weijer, et al. CVPR'07
%   Requires their code from http://cat.uab.es/~joost
%   Just copy mexColorNaming.mex* and w2c.mat

% get image size and obtain patches
[h, w, c] = size(img);
img = double(img);

if w < 10 || h < 10
    error('Image is too small!');
end

patch_coordx = 5:5:w;
patch_coordy = 5:5:h;
[X, Y] = meshgrid(patch_coordy, patch_coordx);
patch = [Y(:), X(:)];
npatches = length(patch);

det = [patch, ones(npatches,1)*2, zeros(npatches,2)];

patch_size = 21;
scalef = 5.0;
sigma = 21 / 4.0;

load 'w2c.mat';
out = mexColorNaming(img(:,:,1), img(:,:,2), img(:,:,3), ...
    det, w2c, patch_size, scalef, sigma);

out = mean(out(:,6:end),1);

end