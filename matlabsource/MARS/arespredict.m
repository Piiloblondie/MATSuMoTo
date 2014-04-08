function Yq = arespredict(model, Xq)
% arespredict
% Predicts output values for the given query points Xq using an ARES model.
%
% Call:
%   Yq = arespredict(model, Xq)
%
% Input:
%   model         : ARES model
%   Xq            : Inputs of query data points (Xq(i,:)), i = 1,...,nq
%
% Output:
%   Yq            : Predicted outputs of the query data points (Yq(i)),
%                   i = 1,...,nq

% =========================================================================
% ARESLab: Adaptive Regression Splines toolbox for Matlab
% Author: Gints Jekabsons (gints.jekabsons@rtu.lv)
% URL: http://www.cs.rtu.lv/jekabsons/
%
% Copyright (C) 2009-2010  Gints Jekabsons
%
% =========================================================================

% Last update: November 9, 2009




if nargin < 2
    error('Too few input arguments.');
end

[mX,nX]=size(Xq); %dimensions of the points where function value should be predicted
[mS,nS]=size(model.minX); %dimesnions of sample sites taken so far 
if nX~=nS %check that both matrices are of the right shape
    Xq=Xq';
end



X = ones(size(Xq,1),length(model.knotdims)+1);
if model.trainParams.cubic
    for i = 1 : length(model.knotdims)
       X(:,i+1) = createbasisfunction(Xq, X, model.knotdims{i}, model.knotsites{i}, ...
                  model.knotdirs{i}, model.parents(i), model.minX, model.maxX, model.t1(i,:), model.t2(i,:));
    end
else
    for i = 1 : length(model.knotdims)
       X(:,i+1) = createbasisfunction(Xq, X, model.knotdims{i}, model.knotsites{i}, ...
                  model.knotdirs{i}, model.parents(i), model.minX, model.maxX);
    end
end
Yq = X * model.coefs;
return
