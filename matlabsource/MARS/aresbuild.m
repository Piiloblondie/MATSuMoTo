function [model, time] = aresbuild(Xtr, Ytr, trainParams, weights, modelOld, verbose)
% aresbuild
% Builds a regression model using the Multivariate Adaptive Regression
% Splines technique.
%
% Call:
%   [model, time] = ...
%             aresbuild(Xtr, Ytr, trainParams, weights, modelOld, verbose)
%
% All the arguments, except the first two, of this function are optional.
% Empty values are also accepted (the corresponding default values will be
% used).
%
% Input:
%   Xtr, Ytr      : Training data cases (Xtr(i,:), Ytr(i)), i = 1,...,n.
%                   Note that it is recommended to pre-scale Xtr values to
%                   [0,1] and to standardize Ytr values.
%   trainParams   : A structure of training parameters for the algorithm.
%                   If not provided, default values will be used (see
%                   function aresparams for details).
%   weights       : A vector of data case weights; if supplied, the
%                   algorithm calculates the sum of squared errors
%                   multiplying the squared residuals by the supplied
%                   weights. The length of weights vector must be the same
%                   as the number of data cases (i.e., n). The weights must
%                   be nonnegative.
%   modelOld      : If here an already built ARES model is provided, no
%                   forward phase will be done. Instead this model will be
%                   taken directly to the backward phase and pruned. This
%                   is useful for fast selection of the "best" penalty
%                   trainParams.c value using Cross-Validation e.g., in
%                   arescvc function.
%   verbose       : Set to false for no verbose. (default value = true)
%
% Output:
%   model         : The built ARES model ? a structure with the following
%                   elements:
%     coefs       : Coefficient vector of the regression model (for the
%                   intercept term and each basis function).
%     knotdims	  : Cell array of indexes of used input variables for each
%                   knot in each basis function.
%     knotsites	  : Cell array of knot sites for each knot and used input
%                   variable in each basis function.
%     knotdirs	  : Cell array of directions (-1 or 1) of the hinge
%                   functions for each used input variable in each basis
%                   function.
%     parents     : Vector of indexes of direct parents for each basis
%                   function (0 if there is no direct parent or it is the
%                   intercept term).
%     trainParams : A structure of training parameters for the algorithm
%                   (the same as in the function input).
%     MSE         : Mean Squared Error of the model in the training data.
%     GCV         : Generalized Cross-Validation (GCV) of the model in the
%                   training data set. The value may also be Inf if model?s
%                   effective number of parameters is larger than or equal
%                   to n.
%     t1          : For piecewise-cubic models only. Matrix of knot sites
%                   for the additional knots on the left of the central
%                   knot.
%     t2          : For piecewise-cubic models only. Matrix of knot sites
%                   for the additional knots on the right of the central
%                   knot.
%     minX        : Vector of minimums for input variables (used for t1 and
%                   t2 placements as well as for model plotting).
%     maxX        : Vector of maximums for input variables (used for t1 and
%                   t2 placements as well as for model plotting).
%     endSpan     : The used value of endSpan.
%   time          : Algorithm execution time (in seconds)

% =========================================================================
% ARESLab: Adaptive Regression Splines toolbox for Matlab
% Author: Gints Jekabsons (gints.jekabsons@rtu.lv)
% URL: http://www.cs.rtu.lv/jekabsons/
%
% Copyright (C) 2009-2010  Gints Jekabsons
%
% =========================================================================

% Please give a reference to the webpage in any publication describing
% research performed using the software e.g., like this:
% Jekabsons G. ARESLab: Adaptive Regression Splines toolbox for Matlab,
% 2010, available at http://www.cs.rtu.lv/jekabsons/

% Last update: April 20, 2010

if nargin < 2
    error('Too few input arguments.');
end

if (isempty(Xtr)) || (isempty(Ytr))
    error('Training data is empty.');
end
if (~isfloat(Xtr)) || (~isfloat(Ytr))
    error('Data type should be floating-point.');
end
[n, d] = size(Xtr); % number of data cases and number of input variables
if size(Ytr,1) ~= n
    error('The number of rows in the matrix and the vector should be equal.');
end
if size(Ytr,2) ~= 1
    error('Ytr should have one column.');
end

if (nargin < 3) || (isempty(trainParams))
    trainParams = aresparams();
end
if trainParams.maxInteractions >= 2
    trainParams_actual_c = trainParams.c;
else
    trainParams_actual_c = 2*trainParams.c/3; % penalty coefficient for additive modelling
end
if (trainParams.cubic) && (trainParams.selfInteractions > 1)
    trainParams.selfInteractions = 1;
    disp('Warning: trainParams.selfInteractions value reverted to 1 due to piecewise-cubic model setting.');
end
if trainParams.cubic
    doCubicFastLevel = trainParams.cubicFastLevel;
    if trainParams.cubicFastLevel > 0
        trainParams.cubic = false; %let's turn it off until the backward phase or the final model
    end
else
    doCubicFastLevel = -1;
end
if trainParams.useMinSpan == 0
    trainParams.useMinSpan = 1; % 1 and 0 is the same here (no endspan)
end
if trainParams.useEndSpan == 0
    trainParams.useEndSpan = 1; % 1 and 0 is the same here (no endspan)
end
if (nargin < 4)
    weights = [];
else
    if (~isempty(weights)) && ...
       ((size(weights,1) ~= n) || (size(weights,2) ~= 1))
        error('weights vector is of wrong size.');
    end
end
wd = diag(weights);
if nargin < 5
    modelOld = [];
end
if (nargin < 6) || (isempty(verbose))
    verbose = true;
end

%if verbose, fprintf('Building ARES model...\n'); end
ws = warning('off');
ttt = tic;

maxIters = floor(trainParams.maxFuncs / 2); % because basis functions are added two at a time
YtrMean = mean(Ytr);
YtrSS = sum((Ytr - YtrMean) .^ 2);
minX = min(Xtr);
maxX = max(Xtr);

if trainParams.useEndSpan < 0
    endSpan = getEndSpan(d); % automatic
else
    endSpan = trainParams.useEndSpan;
end

if isempty(modelOld)
    X = ones(n,1);
    err = 1; % normalized error for the constant model
    model.coefs = YtrMean;
    model.knotdims = {};
    model.knotsites = {};
    model.knotdirs = {};
    model.parents = [];
    model.trainParams = [];
    model.MSE = Inf;
    model.GCV = Inf;
else
    model = modelOld; % modelOld is the initial model
end

if endSpan*2 >= n
    if isempty(modelOld)
        model.MSE = YtrSS / n;
        model.GCV = gcv(model, model.MSE, n, trainParams_actual_c);
        if trainParams.cubic
            model.t1 = [];
            model.t2 = [];
        end
    end
else

    % FORWARD PHASE

    if isempty(modelOld) % no forward phase when modelOld is used

        %if verbose, fprintf('Forward phase  .'); end

        % create sorted lists of data cases for knot placements
        [sortedXtr sortedXtrInd] = sort(Xtr);
        if trainParams.useEndSpan ~= 1
            % throw away data cases at the ends of the intervals
            sortedXtr = sortedXtr(endSpan:end-(endSpan-1),:);
            sortedXtrInd = sortedXtrInd(endSpan:end-(endSpan-1),:);
        end

        if trainParams.cubic
            tmp_t1 = [];
            tmp_t2 = [];
        end
        basisFunctionList = []; % will contain candidate basis functions
        numNewFuncs = 0; % how many basis functions added in the last iteration

        % the main loop of the forward phase
        for depth = 1 : maxIters
            basisFunctionList = createList(basisFunctionList, Xtr, sortedXtr, sortedXtrInd, ...
                                           n, d, model, numNewFuncs, trainParams, endSpan);

            % stop the forward phase if basisFunctionList is empty
            if isempty(basisFunctionList)
                if trainParams.cubic
                    t1 = tmp_t1;
                    t2 = tmp_t2;
                end
                break;
            end

            tmpErr = inf(1,size(basisFunctionList,2));
            tmpCoefs = inf(length(model.coefs)+2, size(basisFunctionList,2));
            Xtmp = zeros(n,size(X,2)+2);
            if ~trainParams.cubic
                Xtmp(:,1:end-2) = X;
            end

            % try all the reflected pairs in the list
            for i = 1 : size(basisFunctionList,2)
                if trainParams.cubic
                    [t1 t2 dif] = findsideknots(model, basisFunctionList{1,i}, basisFunctionList{2,i}, ...
                                  d, minX, maxX, tmp_t1, tmp_t2);
                    Xtmp(:,1:end-2) = X;
                    % update basis functions with the updated side knots
                    for j = 1 : length(model.knotdims)
                        if dif(j)
                            Xtmp(:,j+1) = createbasisfunction(Xtr, Xtmp, model.knotdims{j}, model.knotsites{j}, ...
                                          model.knotdirs{j}, model.parents(j), minX, maxX, t1(j,:), t2(j,:));
                        end
                    end
                    % New basis function
                    dirs = basisFunctionList{3,i};
                    Xtmp(:,end-1) = createbasisfunction(Xtr, Xtmp, basisFunctionList{1,i}, basisFunctionList{2,i}, ...
                                    dirs, basisFunctionList{4,i}, minX, maxX, t1(end,:), t2(end,:));
                    if isnan(Xtmp(1,end-1)), Xtmp(:,end-1) = []; end
                    % Reflected partner
                    dirs(end) = -dirs(end);
                    Xtmp(:,end) = createbasisfunction(Xtr, Xtmp, basisFunctionList{1,i}, basisFunctionList{2,i}, ...
                                  dirs, basisFunctionList{4,i}, minX, maxX, t1(end,:), t2(end,:));
                    if isnan(Xtmp(1,end)), Xtmp(:,end) = []; end
                else
                    % New basis function
                    dirs = basisFunctionList{3,i};
                    Xtmp(:,end-1) = createbasisfunction(Xtr, Xtmp, basisFunctionList{1,i}, ...
                                    basisFunctionList{2,i}, dirs, basisFunctionList{4,i}, minX, maxX);
                    if isnan(Xtmp(1,end-1)), Xtmp(:,end-1) = []; end
                    % Reflected partner
                    dirs(end) = -dirs(end);
                    Xtmp(:,end) = createbasisfunction(Xtr, Xtmp, basisFunctionList{1,i}, ...
                                  basisFunctionList{2,i}, dirs, basisFunctionList{4,i}, minX, maxX);
                    if isnan(Xtmp(1,end)), Xtmp(:,end) = []; end
                end
                [coefs tmpErr(i)] = lreg(Xtmp, Ytr, weights, wd);
                tmpErr(i) = tmpErr(i) / YtrSS;
                if size(Xtmp,2) < size(X,2)+2 % one of the basis functions is not created
                    coefs = [coefs; NaN];
                    Xtmp = [Xtmp zeros(n,1)];
                end
                tmpCoefs(:,i) = coefs;
            end

            [newErr, ind] = min(tmpErr); % find out the best modification

            % stop the forward phase if no correct model was created or
            % if the decrease in error is below the threshold
            if (isnan(newErr)) || (err(end) - newErr < trainParams.threshold)
                if trainParams.cubic
                    t1 = tmp_t1;
                    t2 = tmp_t2;
                end
                break;
            end

            if trainParams.cubic
                [t1 t2 dif] = findsideknots(model, basisFunctionList{1,ind}, basisFunctionList{2,ind}, ...
                              d, minX, maxX, tmp_t1, tmp_t2);
                % update basis functions with the updated side knots
                for j = 1 : length(model.knotdims)
                    if dif(j)
                        X(:,j+1) = createbasisfunction(Xtr, X, model.knotdims{j}, model.knotsites{j}, ...
                                   model.knotdirs{j}, model.parents(j), minX, maxX, t1(j,:), t2(j,:));
                    end
                end
                % Add the new basis function
                dirs = basisFunctionList{3, ind};
                Xn = createbasisfunction(Xtr, X, basisFunctionList{1,ind}, basisFunctionList{2,ind}, ...
                     dirs, basisFunctionList{4,ind}, minX, maxX, t1(end,:), t2(end,:));
                if isnan(Xn(1)), Xn = []; end
                % Add the reflected partner
                dirs(end) = -dirs(end);
                Xn2 = createbasisfunction(Xtr, X, basisFunctionList{1,ind}, basisFunctionList{2,ind}, ...
                      dirs, basisFunctionList{4,ind}, minX, maxX, t1(end,:), t2(end,:));
                if isnan(Xn2(1)), Xn2 = []; end
                X = [X Xn Xn2];
                if ~isempty(Xn) && ~isempty(Xn2) % one of the basis functions is not created
                    t1(end+1,:) = t1(end,:);
                    t2(end+1,:) = t2(end,:);
                end
            else
                dirs = basisFunctionList{3, ind};
                % Add the new basis function
                Xn = createbasisfunction(Xtr, X, basisFunctionList{1,ind}, ...
                     basisFunctionList{2,ind}, dirs, basisFunctionList{4,ind}, minX, maxX);
                if isnan(Xn(1)), Xn = []; end
                % Add the reflected partner
                dirs(end) = -dirs(end);
                Xn2 = createbasisfunction(Xtr, X, basisFunctionList{1,ind}, ...
                      basisFunctionList{2,ind}, dirs, basisFunctionList{4,ind}, minX, maxX);
                if isnan(Xn2(1)), Xn2 = []; end
                X = [X Xn Xn2];
            end

            model.coefs = tmpCoefs(:,ind);

            % add the basis functions to the model
            numNewFuncs = 0;
            dirs = basisFunctionList{3, ind};
            if ~isempty(Xn)
                model.knotdims{end+1,1} = basisFunctionList{1, ind};
                model.knotsites{end+1,1} = basisFunctionList{2, ind};
                model.knotdirs{end+1,1} = dirs;
                model.parents(end+1,1) = basisFunctionList{4, ind};
                numNewFuncs = numNewFuncs + 1;
            else
                model.coefs(end) = [];
            end
            if ~isempty(Xn2)
                dirs(end) = -dirs(end);
                model.knotdims{end+1,1} = basisFunctionList{1, ind};
                model.knotsites{end+1,1} = basisFunctionList{2, ind};
                model.knotdirs{end+1,1} = dirs;
                model.parents(end+1,1) = basisFunctionList{4, ind};
                numNewFuncs = numNewFuncs + 1;
            else
                model.coefs(end) = [];
            end

            %if verbose, fprintf('..'); end

            % stop the forward phase if newErr is too small or if the
            % number of model's coefficients (including the one for the
            % intercept term) in the next iteration is expected to be > n
            err(end+1) = newErr;
            if (newErr < trainParams.threshold) || ...
               (length(model.coefs) + 2 > n)
                break;
            end

            if trainParams.cubic
                tmp_t1 = t1;
                tmp_t2 = t2;
            end
            basisFunctionList(:,ind) = [];
        end % end of the main loop

        %if verbose, fprintf('\n'); end

    end % end of "isempty(modelOld)"

    if isempty(modelOld)
        if (doCubicFastLevel == 1) || ...
           ((doCubicFastLevel >= 2) && (~trainParams.prune))
            % turn the cubic modelling on
            trainParams.cubic = true;
            [t1 t2] = findsideknots(model, [], [], d, minX, maxX, [], []);
            % update all the basis functions
            for i = 1 : length(model.knotdims)
                X(:,i+1) = createbasisfunction(Xtr, X, model.knotdims{i}, model.knotsites{i}, ...
                           model.knotdirs{i}, model.parents(i), minX, maxX, t1(i,:), t2(i,:));
            end
            [model.coefs model.MSE] = lreg(X, Ytr, weights, wd);
            model.MSE = model.MSE / n;
        else
            model.MSE = err(end) * YtrSS / n;
        end
        model.GCV = gcv(model, model.MSE, n, trainParams_actual_c);
        if trainParams.cubic
            model.t1 = t1;
            model.t2 = t2;
        end
    end

    % BACKWARD PHASE

    if trainParams.prune

        %if verbose, fprintf('Backward phase .'); end

        if ~isempty(modelOld) % create basis functions from scratch when modelOld is used
            if (doCubicFastLevel == -1) || (doCubicFastLevel >= 2)
                % create all the basis functions (linear) from scratch
                X = ones(n,length(model.knotdims)+1);
                for i = 1 : length(model.knotdims)
                    X(:,i+1) = createbasisfunction(Xtr, X, model.knotdims{i}, model.knotsites{i}, ...
                               model.knotdirs{i}, model.parents(i), minX, maxX);
                end
                [model.coefs model.MSE] = lreg(X, Ytr, weights, wd);
                model.MSE = model.MSE / n;
                model.GCV = gcv(model, model.MSE, n, trainParams_actual_c);
            else
                % create all the basis functions (cubic) from scratch
                t1 = model.t1;
                t2 = model.t2;
                X = ones(n,length(model.knotdims)+1);
                for i = 1 : length(model.knotdims)
                    X(:,i+1) = createbasisfunction(Xtr, X, model.knotdims{i}, model.knotsites{i}, ...
                               model.knotdirs{i}, model.parents(i), minX, maxX, t1(i,:), t2(i,:));
                end
            end
        end

        models = {model};
        mses = model.MSE;
        gcvs = model.GCV;

        % the main loop of the backward phase
        for j = 1 : length(model.knotdims)
            tmpErr = inf(1, length(model.knotdims));
            tmpCoefs = inf(length(model.coefs)-1, length(model.knotdims));

            % try to delete model's basis functions one at a time
            for k = 1 : length(model.knotdims)
                Xtmp = X;
                Xtmp(:,k+1) = [];
                if trainParams.cubic
                    % create temporary t1, t2, and model with a deleted basis function
                    tmp_t1 = t1;
                    tmp_t1(k,:) = [];
                    tmp_t2 = t2;
                    tmp_t2(k,:) = [];
                    tmp_model.knotdims = model.knotdims;
                    tmp_model.knotdims(k) = [];
                    tmp_model.knotsites = model.knotsites;
                    tmp_model.knotsites(k) = [];
                    tmp_model.knotdirs = model.knotdirs;
                    tmp_model.knotdirs(k) = [];
                    tmp_model.parents = model.parents;
                    tmp_model.parents(k) = [];
                    tmp_model.parents = updateParents(tmp_model.parents, k);
                    [tmp_t1 tmp_t2 dif] = findsideknots(tmp_model, [], [], d, minX, maxX, tmp_t1, tmp_t2);
                    % update basis functions with the updated side knots
                    for i = 1 : length(tmp_model.knotdims)
                        if dif(i)
                            Xtmp(:,i+1) = createbasisfunction(Xtr, Xtmp, tmp_model.knotdims{i}, tmp_model.knotsites{i}, ...
                                          tmp_model.knotdirs{i}, tmp_model.parents(i), minX, maxX, tmp_t1(i,:), tmp_t2(i,:));
                        end
                    end
                end
                [coefs tmpErr(k)] = lreg(Xtmp, Ytr, weights, wd);
                tmpCoefs(:,k) = coefs;
            end

            [dummy, ind] = min(tmpErr); % find out the best modification
            X(:,ind+1) = [];
            model.coefs = tmpCoefs(:,ind);
            model.knotdims(ind) = [];
            model.knotsites(ind) = [];
            model.knotdirs(ind) = [];
            model.parents(ind) = [];
            model.parents = updateParents(model.parents, ind);

            if trainParams.cubic
                t1(ind,:) = [];
                t2(ind,:) = [];
                [t1 t2 dif] = findsideknots(model, [], [], d, minX, maxX, t1, t2);
                % update basis functions with the updated side knots
                for i = 1 : length(model.knotdims)
                    if dif(i)
                        X(:,i+1) = createbasisfunction(Xtr, X, model.knotdims{i}, model.knotsites{i}, ...
                                   model.knotdirs{i}, model.parents(i), minX, maxX, t1(i,:), t2(i,:));
                    end
                end
                model.t1 = t1;
                model.t2 = t2;
            end

            models{end+1} = model;
            mses(end+1) = tmpErr(ind) / n;
            gcvs(end+1) = gcv(model, mses(end), n, trainParams_actual_c);

            %if verbose, fprintf('.'); end
        end % end of the main loop

        if trainParams.maxFinalFuncs >= length(models{1}.coefs)
            [g, ind] = min(gcvs);
        elseif trainParams.maxFinalFuncs > 1
            [g, ind] = min(gcvs(end-trainParams.maxFinalFuncs+1:end));
            ind = ind + length(gcvs) - trainParams.maxFinalFuncs;
        else
            g = gcvs(end);
            ind = length(gcvs);
        end
        model = models{ind};

        if doCubicFastLevel >= 2
            % turn the cubic modelling on
            trainParams.cubic = true;
            [t1 t2] = findsideknots(model, [], [], d, minX, maxX, [], []);
            % update all the basis functions
            X = ones(n,length(model.coefs));
            for i = 1 : length(model.knotdims)
                X(:,i+1) = createbasisfunction(Xtr, X, model.knotdims{i}, model.knotsites{i}, ...
                           model.knotdirs{i}, model.parents(i), minX, maxX, t1(i,:), t2(i,:));
            end
            model.t1 = t1;
            model.t2 = t2;
            [model.coefs model.MSE] = lreg(X, Ytr, weights, wd);
            model.MSE = model.MSE / n;
            model.GCV = gcv(model, model.MSE, n, trainParams_actual_c);
        else
            model.MSE = mses(ind);
            model.GCV = g;
        end

        %if verbose, fprintf('\n'); end

    end % end of "trainParams.prune"

end % end of "if endSpan*2 >= n"

model.trainParams = trainParams;
model.minX = minX;
model.maxX = maxX;
model.endSpan = endSpan;

time = toc(ttt);
if verbose
    %fprintf('Number of basis functions in the final model: %d\n', length(model.coefs));
    maxDeg = 0;
    if length(model.knotdims) > 0
        for i = 1 : length(model.knotdims)
            if length(model.knotdims{i}) > maxDeg
                maxDeg = length(model.knotdims{i});
            end
        end
    end
    %fprintf('Highest degree of interactions in the final model: %d\n', maxDeg);
    %fprintf('Execution time: %0.2f seconds\n', time);
end
warning(ws);
return


%=========================  Auxiliary functions  ==========================

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

function parents = updateParents(parents, deletedInd)
% Updates direct parent indexes after deletion of a basis function.
parents(parents == deletedInd) = 0;
tmp = parents > deletedInd;
parents(tmp) = parents(tmp) - 1;
return

function basisFunctionList = createList(basisFunctionList_old, Xtr, ...
sortedXtr, sortedXtrInd, n, d, model, numNewFuncs, trainParams, endSpan)
% Takes the old list of basis functions and adds new ones according to the
% current model. If the old list is empty, adds only linear basis
% functions. If it is non-empty, adds only basis functions with
% interactions which result from the last numNewFuncs basis functions
% (typically 2, sometimes 1).
% While creating new basis functions, minSpan, endSpan, maxInteractions,
% and selfInteractions are taken into account. Additionally for knot
% placement (and for calculation of minSpan) only those potential knot
% sites are considered for which the parent basis function is larger than
% zero.
% Input argument Xtr should contain all the training data cases; sortedXtr
% should contain all the Xtr cases (except first and last endSpan-1 points)
% sorted in each column separately; sortedXtrInd should contain indexes of
% the sorted cases.

% Create linear basis functions

if (isempty(basisFunctionList_old)) && (numNewFuncs == 0)
    basisFunctionList = cell(3,0);
    counter = 1;
    if trainParams.useMinSpan ~= 1
        % get the list of knot sites allowed due to minSpan
        allowedInd = mod(1:n, getMinSpan(d, n, trainParams.useMinSpan)) == 0;
    end
    for di = 1 : d %for each dimension
        if trainParams.useMinSpan ~= 1
            allowedSites = sortedXtr(:,di);
        else
            allowedSites = unique(sortedXtr(:,di));
        end
        lastknot = Inf;
        % add new unique basis functions to the list (either all or except
        % those which do not fall on the allowed knot sites)
        for i = 1 : size(allowedSites,1)
            if (trainParams.useMinSpan == 1) || ...
               ((allowedInd(i + endSpan-1)) && (lastknot ~= allowedSites(i)))
                lastknot = allowedSites(i);
                basisFunctionList{1, counter} = di;
                basisFunctionList{2, counter} = lastknot;
                basisFunctionList{3, counter} = 1;
                basisFunctionList{4, counter} = 0;
                counter = counter + 1;
            end
        end
    end
    return
end

if (trainParams.maxInteractions < 2) || (numNewFuncs < 1)
    basisFunctionList = basisFunctionList_old;
    return
end

% Create basis functions with interactions

basisFunctionList = basisFunctionList_old;
counter = size(basisFunctionList_old,2) + 1;
start = length(model.knotdims)-(numNewFuncs-1);

% loop through all the basis functions already in the model
for j = start : length(model.knotdims)
    if length(model.knotdims{j}) < trainParams.maxInteractions
        alloweddims = 1 : d;
        if trainParams.selfInteractions <= 1
            % will not consider the already used dimensions
            alloweddims = setdiff(alloweddims, model.knotdims{j});
        else
            for i = 1 : d
                if length(find(model.knotdims{j} == i)) >= trainParams.selfInteractions
                    alloweddims = setdiff(alloweddims, i);
                end
            end
        end
        if isempty(alloweddims)
            continue
        end

        if trainParams.useMinSpan ~= 1

            nonzero = listNonZero(Xtr, model.knotdims{j}, model.knotsites{j}, model.knotdirs{j});
            minSpan = getMinSpan(d, length(find(nonzero)), trainParams.useMinSpan);
            if ~isfinite(minSpan)
                continue
            end
            allowedInd = mod(1:n, minSpan) == 0;
            for di = alloweddims %for each dimension
                lastknot = Inf;
                % add new unique basis functions to the list (all except
                % those which do not fall on the allowed knot sites and
                % except those for which the parent basis function is zero
                % on the knot site)
                ind = sortedXtrInd(nonzero(sortedXtrInd(:,di)),di);
                ind = ind(allowedInd(ind))';
                for i = ind
                    if lastknot ~= Xtr(i,di)
                        lastknot = Xtr(i,di);
                        basisFunctionList{1, counter} = [model.knotdims{j} di];
                        basisFunctionList{2, counter} = [model.knotsites{j} lastknot];
                        basisFunctionList{3, counter} = [model.knotdirs{j} 1];
                        basisFunctionList{4, counter} = j;
                        counter = counter + 1;
                    end
                end
            end

        else

            nonzero = listNonZero(Xtr, model.knotdims{j}, model.knotsites{j}, model.knotdirs{j});
            for di = alloweddims %for each dimension
                % add new unique basis functions to the list (all except
                % those for which the parent basis function is zero on
                % the knot site)
                allowedSites = unique(sortedXtr(nonzero(sortedXtrInd(:,di)),di));
                for i = 1 : size(allowedSites,1)
                    basisFunctionList{1, counter} = [model.knotdims{j} di];
                    basisFunctionList{2, counter} = [model.knotsites{j} allowedSites(i)];
                    basisFunctionList{3, counter} = [model.knotdirs{j} 1];
                    basisFunctionList{4, counter} = j;
                    counter = counter + 1;
                end
            end

        end

    end
end
return

function nonzero = listNonZero(Xtr, knotdims, knotsites, knotdirs)
% Lists nonzero (according to the parent basis function) sites where knots
% may be placed.
nonzero = true(size(Xtr,1),1);
for j = 1 : size(Xtr,1)
    for i = 1 : length(knotdims)
        z = Xtr(j,knotdims(i)) - knotsites(i);
        if ((z >= 0) && (knotdirs(i) < 0)) || ...
           ((z <= 0) && (knotdirs(i) > 0))
            nonzero(j) = false;
            break;
        end
    end
end
return

function s = getEndSpan(d)
% Calculation of endSpan so that potential knot sites that are too close to
% the ends of data intervals are not considered.
%s = floor(3 - log2(0.05/d));
s = floor(7.32193 + log(d) / 0.69315); % precomputed version
if s < 1, s = 1; end
return

function s = getMinSpan(d, nz, param)
% Calculation of minSpan so that only those potential knot sites are
% considered which are at least minSpan apart. This increases resistance to
% runs of correlated noise.
% nz is the number of potential knot sites where the parent basis function
% is nonzero.
if nz == 0
    s = Inf;
else
    if param < 0 % automatic
        %s = floor(-log2(-log(1-0.05)/(d*nz)) / 2.5);
        s = floor((2.9702 + log(d*nz)) / 1.7329); % precomputed version
    else
        s = param;
        if s > nz
            s = Inf;
        end
    end
    if s < 1, s = 1; end
end
return

function [coefs err] = lreg(x, y, w, wd)
% Linear regression (unweighted and weighted)
if isempty(wd)
    coefs = (x' * x) \ (x' * y);
    err = sum((y-x*coefs).^2);
else
    x_wd = x' * wd;
    coefs = (x_wd * x) \ (x_wd * y);
    err = sum((y-x*coefs).^2.*w);
end
return
