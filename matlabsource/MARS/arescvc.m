function cBest = arescvc(X, Y, trainParams, cTry, weights, k, shuffle, verbose)
% arescvc
% Finds the "best" value for penalty c from a set of candidate values using
% k-fold Cross-Validation and MSE.
%
% Call:
%   cBest = arescv(X, Y, trainParams, cTry, weights, k, shuffle, verbose)
%
% All the arguments, except the first three, of this function are optional.
% Empty values are also accepted (the corresponding default values will be
% used).
%
% Input:
%   X, Y          : Data cases (X(i,:), Y(i)), i = 1,...,n
%   trainParams   : See function aresbuild.
%   cTry          : A set of candidate values for c. (default = 1:5)
%   weights       : See function aresbuild.
%   k             : Value of k for k-fold Cross-Validation. The typical
%                   values are 5 or 10. For Leave-One-Out Cross-Validation
%                   set k equal to n. (default value = 10)
%   shuffle       : Whether to shuffle the order of the data cases before
%                   performing Cross-Validation. Note that the random seed
%                   value can be controlled externally before calling
%                   arescvc. (default value = true)
%   verbose       : Set to false for no verbose. (default value = true)
%
% Output:
%   cBest         : The "best" value for penalty c.

% =========================================================================
% ARESLab: Adaptive Regression Splines toolbox for Matlab
% Author: Gints Jekabsons (gints.jekabsons@rtu.lv)
% URL: http://www.cs.rtu.lv/jekabsons/
%
% Copyright (C) 2009-2010  Gints Jekabsons
%
% =========================================================================

% Last update: April 20, 2010

if nargin < 2
    error('Too few input arguments.');
end

if (isempty(X)) || (isempty(Y))
    error('Data is empty.');
end
[n, d] = size(X); % number of data cases and number of input variables
if size(Y,1) ~= n
    error('The number of rows in the matrix and the vector should be equal.');
end
if size(Y,2) ~= 1
    error('The vector Y should have one column.');
end

if (nargin < 3) || (isempty(trainParams))
    trainParams = aresparams();
else
    if ~trainParams.prune
        error('Model pruning is disabled (does not make sense because penalty c is used only in pruning).');
    end
end

if (nargin < 4) || (isempty(cTry))
    cTry = 1:5;
end
cNum = length(cTry);
if (nargin < 5)
    weights = [];
end
if (nargin < 6) || (isempty(k))
    k = 10;
end
if k < 2
    error('k should not be smaller than 2.');
end
if k > n
    error('k should not be larger than the number of data cases.');
end
if (nargin < 7) || (isempty(shuffle))
    shuffle = true;
end
if (nargin < 8) || (isempty(verbose))
    verbose = true;
end

if shuffle
    ind = randperm(n); % shuffle the data
else
    ind = 1 : n;
end

% divide the data into k subsets
minsize = floor(n / k);
remainder = n - minsize * k;
sizes = zeros(k, 1);
for i = 1 : k
    sizes(i) = minsize;
    if remainder > 0
        sizes(i) = sizes(i) + 1;
        remainder = remainder - 1;
    end
end
offsets = ones(k, 1);
for i = 2 : k
    offsets(i) = offsets(i-1) + sizes(i-1);
end

% perform the training and testing
MSE = zeros(cNum,1);
for i = 1 : k
    Xtr = zeros(n-sizes(k-i+1), d);
    Ytr = zeros(n-sizes(k-i+1), 1);
    currsize = 0;
    for j = 1 : k
        if k-i+1 ~= j
            Xtr(currsize+1 : currsize+1+sizes(j)-1, :) = X(ind(offsets(j):offsets(j)+sizes(j)-1), :);
            Ytr(currsize+1 : currsize+1+sizes(j)-1, 1) = Y(ind(offsets(j):offsets(j)+sizes(j)-1), 1);
            currsize = currsize + sizes(j);
        end
    end
    Xtst = X(ind(offsets(k-i+1):offsets(k-i+1)+sizes(k-i+1)-1), :);
    Ytst = Y(ind(offsets(k-i+1):offsets(k-i+1)+sizes(k-i+1)-1), 1);
    if verbose
        disp(['Fold #' num2str(i)]);
    end
    trainParams.prune = false;
    modelBase = aresbuild(Xtr, Ytr, trainParams, weights, [], false); % unpruned base model
    trainParams.prune = true;
    for j = 1 : cNum
        trainParams.c = cTry(j);
        model = aresbuild(Xtr, Ytr, trainParams, weights, modelBase, false); % pruned model
        MSE(j) = MSE(j) + arestest(model, Xtst, Ytst);
    end
end

MSE = MSE ./ k;
if verbose
    for j = 1 : cNum
        disp(['c = ' num2str(cTry(j)) '  MSE = ' num2str(MSE(j))]);
    end
end
[dummy ind] = min(MSE);
cBest = cTry(ind);
return
