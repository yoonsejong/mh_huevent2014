function [fcm_model] = fuzzyFuncLearnItten()
%   FUZZY_FUNC_LEARN
%
%   Do Fuzzy C-Means clustering on affective image classification dataset
%   following Machajdik and Hanbury, MM'13.
%

%% Please uncomment following lines to start from segmentation (takes a little time)
% path_to_AIC_dataset = 'D:/.../AffectiveImageClassification/';
% 
% [seg1lab, seg1avg, seg1var] = watershedSegmentBatch([path_to_AIC_dataset 'testImages_abstract'], 'jpg');
% [seg2lab, seg2avg, seg2var] = watershedSegmentBatch([path_to_AIC_dataset 'testImages_artphoto'], 'jpg');
% 
% save('AIC_watershed_segmentation.mat', 'seg1lab', 'seg1avg', ...
%     'seg2lab', 'seg2avg', 'seg1var', 'seg2var');   

%% Please uncomment following lines to re-train Fuzzy membership function
% 
% load('AIC_watershed_segmentation.mat');
% 
% nSeg1 = zeros(size(seg1avg,1),1);
% for idx = 1 : size(nSeg1,1), nSeg1(idx) = size(seg1avg{idx},1); end;
% nSeg2 = zeros(size(seg2avg,1),1);
% for idx = 1 : size(nSeg2,1), nSeg2(idx) = size(seg2avg{idx},1); end;
% 
% seg1avgMat = cell2mat(seg1avg);
% seg2avgMat = cell2mat(seg2avg);
% seg3avgMat = [seg1avgMat; seg2avgMat];
% 
% % (1) exponent m for U
% % (2) max iterations
% % (3) min improvement (tolerence)
% % (4) show iteration info or not
% fcmOptions = [2.0, 100, 1e-5, 1];
% 
% % Do Fuzzy C Means
% [Cs1, Us1, ~] = fcm(seg1avgMat(:,2), 3, fcmOptions); % 3 level saturation (S)
% [Cy1, Uy1, ~] = fcm(seg1avgMat(:,3), 5, fcmOptions); % 5 level luminance (Y)
% Us1 = Us1';
% Uy1 = Uy1';
% 
% [Cs2, Us2, ~] = fcm(seg2avgMat(:,2), 3, fcmOptions); % 3 level saturation (S)
% [Cy2, Uy2, ~] = fcm(seg2avgMat(:,3), 5, fcmOptions); % 5 level luminance (Y)
% Us2 = Us2';
% Uy2 = Uy2';
% 
% [Cs3, Us3, ~] = fcm(seg3avgMat(:,2), 3, fcmOptions); % 3 level saturation (S)
% [Cy3, Uy3, ~] = fcm(seg3avgMat(:,3), 5, fcmOptions); % 5 level luminance (Y)
% Us3 = Us3';
% Uy3 = Uy3';
% 
% % compute size of each segments
% seg1size = cell(size(seg1lab,1),1);
% for idf = 1 : size(seg1lab,1)
%     segs1 = unique(seg1lab{idf});
%     seg1size{idf} = zeros(length(segs1),1);
%     for ids = 1 : length(segs1)
%         seg1size{idf}(ids) = sum(seg1lab{idf} == segs1(ids));
%     end
% end
% seg2size = cell(size(seg2lab,1),1);
% for idf = 1 : size(seg2lab,1)
%     segs2 = unique(seg2lab{idf});
%     seg2size{idf} = zeros(length(segs2),1);
%     for ids = 1 : length(segs2)
%         seg2size{idf}(ids) = sum(seg2lab{idf} == segs2(ids));
%     end
% end
% 
% save('AIC_fcm_itten.mat', ...
%     'seg1size', 'nSeg1', 'Cs1', 'Cy1', 'Cs2', 'Cy2', 'Cs3', 'Cy3', ...
%     'seg2size', 'nSeg2', 'Us1', 'Uy1', 'Us2', 'Uy2', 'Us3', 'Uy3', 'fcmOptions');

fcm_model = load('AIC_fcm_itten.mat');

end
