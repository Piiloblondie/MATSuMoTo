function aresplot(model, minX, maxX, vals, gridSize)
% aresplot
% Plots surface of an ARES model.
%
% Call:
%   aresplot(model, minX, maxX, vals, gridSize)
%
% All the arguments, except the first one, of this function are optional.
% Empty values are also accepted (the corresponding default values will be
% used).
%
% Input:
%   model         : ARES model
%   minX, maxX    : User defined minimum and maximum values for each input
%                   variable (this is the same type of data as in
%                   model.minX and model.maxX). If not supplied, the
%                   model.minX and model.maxX values will be used.
%   vals          : Only used when the number of input variables is larger
%                   than 2. This is a vector of fixed values for all the
%                   input variables except the two varied in the plot. The
%                   two varied variables are identified in vals using NaN
%                   values. By default the two first variables will be
%                   varied and all the other will be fixed at (max-min)/2.
%   gridSize      : Grid size. (default value = 50)

% =========================================================================
% ARESLab: Adaptive Regression Splines toolbox for Matlab
% Author: Gints Jekabsons (gints.jekabsons@rtu.lv)
% URL: http://www.cs.rtu.lv/jekabsons/
%
% Copyright (C) 2009-2010  Gints Jekabsons
%
% =========================================================================

% Last update: November 9, 2009

if nargin < 1
    error('Too few input arguments.');
end
if (nargin < 2) || (isempty(minX))
    minX = model.minX;
else
    if length(minX) ~= length(model.minX)
        error('Vector minX is of wrong size.');
    end
end
if (nargin < 3) || (isempty(maxX))
    maxX = model.maxX;
else
    if length(maxX) ~= length(model.maxX)
        error('Vector maxX is of wrong size.');
    end
end
if (nargin < 5) || (isempty(gridSize))
    gridSize = 50;
end

if length(model.minX) < 2
    step = (maxX - minX) / gridSize;
    X = [minX:step:maxX]';
    plot(X, arespredict(model, X));
    return
end

if (nargin < 4) || (isempty(vals))
    ind1 = 1; vals(1) = NaN;
    ind2 = 2; vals(2) = NaN;
    if length(minX) > 2
        for i = 3 : length(minX)
            vals(i) = (maxX(i) - minX(i)) / 2;
        end
    end
else
    if length(minX) ~= length(vals)
        error('Vector vals is of wrong size.');
    end
    tmp = 0;
    for i = 1 : length(vals)
        if isnan(vals(i))
            if tmp == 0
                ind1 = i;
                tmp = 1;
            elseif tmp == 1
                ind2 = i;
                tmp = 2;
            else
                tmp = 3;
                break;
            end
        end
    end
    if tmp ~= 2
        error('Vector vals should contain NaN exactly two times.');
    end
end

step1 = (maxX(ind1) - minX(ind1)) / gridSize;
step2 = (maxX(ind2) - minX(ind2)) / gridSize;

[X1,X2] = meshgrid(minX(ind1):step1:maxX(ind1), minX(ind2):step2:maxX(ind2));
XX1 = reshape(X1, numel(X1), 1);
XX2 = reshape(X2, numel(X2), 1);

X = zeros(size(XX1,1), length(minX));
X(:,ind1) = XX1;
X(:,ind2) = XX2;
for i = 1 : length(minX)
    if (i ~= ind1) && (i ~= ind2)
        X(:,i) = vals(i);
    end
end

YY = arespredict(model, X);
Y = reshape(YY, size(X1,1), size(X2,2));
surfc(X1, X2, Y);
return
