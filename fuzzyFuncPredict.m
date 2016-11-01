function [Us, Uy, As, Ay] = fuzzyFuncPredict(model, avgImgSegs)
%   FUZZY_FUNC_PREDICT
%
%   Do one fuzzy C-means iteration WITHOUT updating the model centers.
%   Thus we only predict the membership function for the query samples.
%
% INPUT
%   model      : Trained Fuzzy C-Means centers with parameters
%   avgImgSegs : NSegs x 3 matrix containing mean HSY colors of segments or
%                you can directly specify a vectorized HSY image
%
% OUTPUT
%   Us, Uy     : NSeg x K matrices of predicted membership function
%   As, Ay     : NSeg x 1 vectors of predicted membership of each sample

if isfield(model, 'Cs3')
    nCs = length(model.Cs3);
    nCy = length(model.Cy3);
    
    Cs = model.Cs3;
    Cy = model.Cy3;
else
    nCs = length(model.Cs);
    nCy = length(model.Cy);
    
    Cs = model.Cs;
    Cy = model.Cy;
end
m = model.fcmOptions(1);

dist = distfcm(Cs, avgImgSegs(:,2)); % this will return nCs x NSegs
tmp = dist .^ (-2 / (m-1));
Us = tmp ./ (ones(nCs, 1) * sum(tmp));

dist = distfcm(Cy, avgImgSegs(:,3));
tmp = dist .^ (-2 / (m-1));
Uy = tmp ./ (ones(nCy, 1) * sum(tmp));

Us = Us';
Uy = Uy';

if nargout > 2
    [~, As] = max(Us, [], 2);
    [~, Ay] = max(Uy, [], 2);
end

end