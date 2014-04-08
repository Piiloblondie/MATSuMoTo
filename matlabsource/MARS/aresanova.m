function aresanova(model, Xtr, Ytr)
% aresanova
% Performs ANOVA decomposition (see Sections 3.5 and 4.3 of the original
% paper of Jerome Friedman (Friedman 1991) for details) of the given ARES
% model and reports the results.
%
% Call:
%   aresanova(model, Xtr, Ytr)
%
% Input:
%   model         : ARES model
%   Xtr, Ytr      : Training data cases (Xtr(i,:), Ytr(i)), i = 1,...,n.

% =========================================================================
% ARESLab: Adaptive Regression Splines toolbox for Matlab/Octave
% Author: Gints Jekabsons (gints.jekabsons@rtu.lv)
% URL: http://www.cs.rtu.lv/jekabsons/
%
% Copyright (C) 2009-2010  Gints Jekabsons
%
% =========================================================================

% Last update: May 28, 2010

if nargin ~= 3
    error('The number of input arguments should be axactly three.');
end
if (isempty(Xtr)) || (isempty(Ytr))
    error('Data is empty.');
end
n = size(Xtr,1);
if size(Ytr,1) ~= n
    error('The number of rows in the matrix and the vector should be equal.');
end
if size(Ytr,2) ~= 1
    error('Ytr should have one column.');
end

if model.trainParams.cubic
    fprintf('Type: piecewise-cubic\n');
else
    fprintf('Type: piecewise-linear\n');
end
fprintf('GCV: %0.3f\n', model.GCV);
fprintf('Total number of basis functions: %d\n', length(model.coefs));
fprintf('Total effective number of parameters: %0.1f\n', ...
        length(model.coefs) + model.trainParams.c * length(model.knotdims) / 2);
fprintf('ANOVA decomposition:\n');
fprintf('Func.\t\tSTD\t\t\tGCV\t\t#basis\t#params\t\tvariable(s)\n');
counter = 0;
for i = 1 : model.trainParams.maxInteractions
    combs = nchoosek(1:length(model.minX),i);
    for j = 1 : size(combs,1)
        [modelReduced usedBasis] = aresanovareduce(model, combs(j,:), true);
        if length(usedBasis) > 0
            counter = counter + 1;
            fprintf('%d\t\t', counter);
            fprintf('%7.3f\t', std(arespredict(modelReduced, Xtr))); % standard deviation of the ANOVA function
            modelReduced = deleteBasis(model, usedBasis);
            Yq = arespredict(modelReduced, Xtr);
            MSE = mean((Ytr - Yq).^2);
            fprintf('%11.3f\t\t', gcv(modelReduced, MSE, n, model.trainParams.c)); % GCV when the basis functions are deleted
            fprintf('%6d\t', length(usedBasis)); % the number of basis functions for that ANOVA function
            fprintf('%7.1f\t\t', length(usedBasis) + model.trainParams.c * length(usedBasis) / 2); % effective parameters
            fprintf('%d ', combs(j,:)); % used variables
            fprintf('\n');
        end
    end
end
return

function modelReduced = deleteBasis(model, nums)
modelReduced = model;
for i = length(nums) : -1 : 1
    modelReduced.coefs(nums(i)+1) = [];
    modelReduced.knotdims(nums(i)) = [];
    modelReduced.knotsites(nums(i)) = [];
    modelReduced.knotdirs(nums(i)) = [];
    modelReduced.parents(nums(i)) = [];
    modelReduced.t1(nums(i),:) = [];
    modelReduced.t2(nums(i),:) = [];
end
modelReduced.parents(:) = 0;
if modelReduced.trainParams.cubic
    % correct the side knots for the modified cubic model
    [modelReduced.t1 modelReduced.t2] = ...
    findsideknots(modelReduced, [], [], size(modelReduced.t1,2), modelReduced.minX, modelReduced.maxX, [], []);
end
return

function g = gcv(model, MSE, n, c)
% Calculates GCV from model complexity, its Mean Squared Error, number of
% data cases n, and penalty coefficient c.
enp = length(model.coefs) + c * length(model.knotdims) / 2; % model's effective number of parameters
if enp >= n
    g = Inf;
else
    p = 1 - enp / n;
    g = MSE / (p * p);
end
return
