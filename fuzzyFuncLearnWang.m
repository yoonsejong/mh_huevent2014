function fcm_model = fuzzyFuncLearnWang()
%   FUZZY_FUNC_LEARN_WANG
%
%   Do Fuzzy C-Means clustering on affective image classification dataset
%   following W. Wang, et al. IEEE Intl. Conf. Sys, Man, Cyber. 2006.
%

%% Generate HSY image from AIC dataset. Creates an 1.8 GB file!
% path_to_AIC_dataset = 'D:/.../AffectiveImageClassification/';
% 
% files1 = dir([path_to_AIC_dataset 'testImages_abstract/*.jpg']);
% files2 = dir([path_to_AIC_dataset 'testImages_artphoto/*.jpg']);
% 
% nFiles1 = length(files1);
% nFiles2 = length(files2);
% 
% bigH = cell(nFiles1+nFiles2,1);
% bigS = cell(nFiles1+nFiles2,1);
% bigY = cell(nFiles1+nFiles2,1);
% for idf = 1:nFiles1
%     img = imread([path_to_AIC_dataset 'testImages_abstract/' files1(idf).name]);
%     if ismatrix(img)
%         img = repmat(img, [1, 1, 3]);
%     end
%     img = rgb2hsy(img);
%     [h, w, c] = size(img);
%     img = reshape(img, [h*w, c]);
%     bigH{idf} = img(:,1);
%     bigS{idf} = img(:,2);
%     bigY{idf} = img(:,3);
%     fprintf('Processed %d / %d files.\n', idf, nFiles1+nFiles2);
% end
% for idf = 1+nFiles1:nFiles2+nFiles1
%     img = imread([path_to_AIC_dataset 'testImages_artphoto/' files2(idf-nFiles1).name]);
%     if ismatrix(img)
%         img = repmat(img, [1, 1, 3]);
%     end
%     img = rgb2hsy(img);
%     [h, w, c] = size(img);
%     img = reshape(img, [h*w, c]);
%     bigH{idf} = img(:,1);
%     bigS{idf} = img(:,2);
%     bigY{idf} = img(:,3);
%     fprintf('Processed %d / %d files.\n', idf, nFiles1+nFiles2);
% end
% bigH = cell2mat(bigH);
% bigS = cell2mat(bigS);
% bigY = cell2mat(bigY);
% 
% save('AIC_convert_hsy.mat', 'bigH', 'bigS', 'bigY');   

%% Please uncomment following lines to re-train Fuzzy membership function
% 
% load('AIC_convert_hsy.mat');
% 
% % We randomly sample 3 million pixels from the dataset
% clearvars bigH;
% randIdx = randperm(size(bigS,1));
% randIdx = randIdx(1:3000000);
% bigS = sort(bigS);
% bigY = sort(bigY);
% bigS = bigS(randIdx);
% bigY = bigY(randIdx);
% 
% % (1) exponent m for U
% % (2) max iterations
% % (3) min improvement (tolerence)
% % (4) show iteration info or not
% fcmOptions = [2.0, 100, 1e-5, 1];
% 
% % Do Fuzzy C Means
% [Cs, ~, ~] = fcm(bigS, 2, fcmOptions); % 2 level saturation (S)
% [Cy, ~, ~] = fcm(bigY, 5, fcmOptions); % 5 level luminance (Y)
% 
% Cs = sort(Cs);
% Cy = sort(Cy);
% 
% save('AIC_fcm_wang.mat', 'Cs', 'Cy', 'fcmOptions');

fcm_model = load('AIC_fcm_wang.mat');

end

function [c] = fuzzyFuncLearnWang_impl(x)
% This is an implementation of the algorithm described in Want, et al.
% We merely use Fuzzy C-Means instead.

n = length(x);

% init c0 = min{x1, x2, ..., xn}
c0 = min(x);

% init c6 = max{x1, x2, ..., xn}
c6 = max(x);

% init c1 ... c5
c = zeros(5,1);
for j = 1 : 5
    c(j) = c0 + j * (c6 + c0) / 6;
end

iter = 0;
max_iter = 10000;
tol = 1e-5;
relerr = inf;
oldc = c;
while relerr >= tol && iter < max_iter
    % init and calculate U
    U = zeros(n, 5);

    for i = 1 : n
        for j = 1 : 5
            if x(i) <= c(1)
                U(i,1) = 1; 
                U(i,2:end) = 0; 
            end
            if x(i) > c(5)
                U(i,1:(end-1)) = 0;
                U(i,5) = 1;
            end
            
            if j == 5
                c_j_plus_one = c6;
            else
                c_j_plus_one = c(j+1);
            end
            
            if c(j) < x(i) && x(i) <= c_j_plus_one
                U(i,j) = (c_j_plus_one - x(i)) / (c_j_plus_one - c(j));
                U(i,j+1) = 1 - U(i,j);
                U(i,1:(j-1)) = 1 - U(i,j);
                U(i,(j+1):end) = 1 - U(i,j);
            end
        end
    end

    % update c
    for j = 1 : 5
        c(j) = (U(:,j)' * x) / sum(U(:,j));
    end

    % compute relative error
    relerr = norm(c - oldc) / norm(c);
    
    % increase counter
    iter = iter + 1;
end

end
