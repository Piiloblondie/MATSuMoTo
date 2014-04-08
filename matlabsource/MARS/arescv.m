function [avgMSE, avgRMSE, avgRRMSE, avgR2, avgTime] = arescv(X, Y, ...
trainParams, weights, k, shuffle, cvc_cTry, cvc_k, verbose)
% arescv
% Tests ARES performance using k-fold Cross-Validation.
%
% Call:
%   [avgMSE, avgRMSE, avgRRMSE, avgR2, avgTime] = arescv(X, Y, ...
%              trainParams, weights, k, shuffle, cvc_cTry, cvc_k, verbose)
%
% All the arguments, except the first two, of this function are optional.
% Empty values are also accepted (the corresponding default values will be
% used).
%
% Input:
%   X, Y          : Data cases (X(i,:), Y(i)), i = 1,...,n
%   trainParams   : See function aresbuild.
%   weights       : See function aresbuild.
%   k             : Value of k for k-fold Cross-Validation. The typical
%                   values are 5 or 10. For Leave-One-Out Cross-Validation
%                   set k equal to n. (default value = 10)
%   shuffle       : Whether to shuffle the order of the data cases before
%                   performing Cross-Validation. Note that the random seed
%                   value can be controlled externally before calling
%                   arescv. (default value = true)
%   cvc_cTry, cvc_k : cTry and k values for arescvc function. Supply these
%                   values if you want to perform another Cross-Validation
%                   for finding the "best" penalty c value in each
%                   iteration of the outer Cross-Validation loop of arescv.
%                   (default values = [], meaning that a fixed c is used)
%   verbose       : Set to false for no verbose. (default value = true)
%
% Output:
%   avgMSE        : Average Mean Squared Error
%   avgRMSE       : Average Root Mean Squared Error
%   avgRRMSE      : Average Relative Root Mean Squared Error
%   avgR2         : Average Coefficient of Determination
%   avgTime       : Average execution time

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
end
if (nargin < 4)
    weights = [];
end
if (nargin < 5) || (isempty(k))
    k = 10;
end
if k < 2
    error('k should not be smaller than 2.');
end
if k > n
    error('k should not be larger than the number of data cases.');
end
if (nargin < 6) || (isempty(shuffle))
    shuffle = true;
end
if (nargin < 7) || (isempty(cvc_cTry))
    cvc_cTry = [];
end
if (nargin < 8) || (isempty(cvc_k))
    cvc_k = [];
end
if (nargin < 9) || (isempty(verbose))
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

% perform the training and testing k times
MSE = NaN(k,1);
RMSE = NaN(k,1);
RRMSE = NaN(k,1);
R2 = NaN(k,1);
time = NaN(k,1);
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

    if isempty(cvc_k)
        % perform just model building
        [model, time(i)] = aresbuild(Xtr, Ytr, trainParams, weights, [], verbose);
    else
        % perform finding the "best" value for trainParams.c and then build
        % a model using the "best" value
        tic;
        trainParams.c = arescvc(Xtr, Ytr, trainParams, cvc_cTry, weights, cvc_k, false, false);
        model = aresbuild(Xtr, Ytr, trainParams, weights, [], verbose);
        time(i) = toc;
    end
    [MSE(i), RMSE(i), RRMSE(i), R2(i)] = arestest(model, Xtst, Ytst);
end

avgMSE = mean(MSE);
avgRMSE = mean(RMSE);
avgRRMSE = mean(RRMSE);
avgR2 = mean(R2);
avgTime = mean(time);
return
