function eq = areseq(model, precision)
% areseq
% Outputs the ARES model in an explicit mathematical form (useful e.g., for
% deployment of built ARES models in other software).
%
% Call:
%   eq = areseq(model, precision)
%   eq = areseq(model)
%
% Input:
%   model         : ARES model
%   precision     : Number of digits in the model coefficients and knot
%                   sites.
%
% Output:
%   eq            : A cell array of equations for individual basis
%                   functions and the main model. Note that the outputted
%                   equations will have piecewise-linear form even if the
%                   model itself is really piecewise-cubic. The model
%                   coefficients however will be those from the original
%                   piecewise-cubic model.

% =========================================================================
% ARESLab: Adaptive Regression Splines toolbox for Matlab
% Author: Gints Jekabsons (gints.jekabsons@rtu.lv)
% URL: http://www.cs.rtu.lv/jekabsons/
%
% Copyright (C) 2009-2010  Gints Jekabsons
%
% =========================================================================

% Last update: April 9, 2010

if nargin < 1
    error('Too few input arguments.');
end
if (nargin < 2) || (isempty(precision))
    precision = 15;
end
if model.trainParams.cubic
    disp('Warning: The model is piecewise-cubic but the basis functions will be shown as piecewise-linear.');
end

p = ['%.' num2str(precision) 'g'];
eq = {};

% compose the individual basis functions
for i = 1 : length(model.knotdims)
    tmp = ['BF' num2str(i) ' ='];
    if model.parents(i) > 0
        tmp = [tmp ' BF' num2str(model.parents(i)) ' *'];
        start = length(model.knotdims{i});
    else
        start = 1;
    end
    for j = start : length(model.knotdims{i})
        if model.knotdirs{i}(j) > 0
            if model.knotsites{i}(j) >= 0
                m = '-';
            else
                m = '';
            end
            tmp = [tmp ' max(0, x' num2str(model.knotdims{i}(j),p) ' ' m num2str(model.knotsites{i}(j),p) ')'];
        else
            tmp = [tmp ' max(0, ' num2str(model.knotsites{i}(j),p) ' -x' num2str(model.knotdims{i}(j),p) ')'];
        end
        if j < length(model.knotdims{i})
            tmp = [tmp ' *'];
        end
    end
    disp(tmp);
    eq{end+1} = tmp;
end

% compose the summation
tmp = ['y = ' num2str(model.coefs(1),p)];
for i = 1 : length(model.knotdims)
    if model.coefs(i+1) >= 0
        tmp = [tmp ' +'];
    else
        tmp = [tmp ' '];
    end
    tmp = [tmp num2str(model.coefs(i+1),p) '*BF' num2str(i)];
end
disp(tmp);
eq{end+1} = tmp;
return
