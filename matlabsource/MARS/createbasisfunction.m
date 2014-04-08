function Xn = createbasisfunction(X, Xtmp, knotdims, knotsites, knotdirs, parent, minX, maxX, t1, t2)
% Creates a list of response values of a defined basis function for either
% a piecewise-linear or piecewise-cubic model.

% =========================================================================
% ARESLab: Adaptive Regression Splines toolbox for Matlab
% Author: Gints Jekabsons (gints.jekabsons@rtu.lv)
% URL: http://www.cs.rtu.lv/jekabsons/
%
% Copyright (C) 2009-2010  Gints Jekabsons
%
% =========================================================================

% Last update: November 9, 2009

n = size(X,1);

% if the basis function will always be zero, return NaN
if ((knotdirs(end) > 0) && (knotsites(end) >= maxX(knotdims(end)))) || ...
   ((knotdirs(end) < 0) && (knotsites(end) <= minX(knotdims(end))))
    Xn = NaN(n,1);
    return
end

% if the basis function has a direct parent in the model, use it to speed-up the calculations
if parent > 0
    Xn = Xtmp(:,parent+1);
    start = length(knotdims);
else
    Xn = ones(n,1);
    start = 1;
end

if nargin < 10 % piecewise-linear

    for i = start : length(knotdims)
        if knotdirs(i) > 0
            z = X(:,knotdims(i)) - knotsites(i);
        else
            z = knotsites(i) - X(:,knotdims(i));
        end
        Xn = Xn .* max(0,z);
    end

else % piecewise-cubic

    Xx = zeros(n,1);
    for i = start : length(knotdims)
        % if the knot is on the very edge, treat the basis function as linear
        if (knotdirs(i) > 0) && (knotsites(i) <= minX(knotdims(i)))
            Xx(:) = X(:,knotdims(i)) - knotsites(i);
        elseif (knotdirs(i) < 0) && (knotsites(i) >= maxX(knotdims(i)))
            Xx(:) = knotsites(i) - X(:,knotdims(i));
        else
            tt1 = t1(knotdims(i));
            tt2 = t2(knotdims(i));
            if knotdirs(i) > 0
                % p = (2*tt2 + tt1 - 3*knotsites(i)) / (tt2 - tt1)^2; % explicit version
                % r = (2*knotsites(i) - tt2 - tt1) / (tt2 - tt1)^3; % explicit version
                tt2_tt1 = tt2 - tt1;
                tt2_tt1_sq = tt2_tt1 * tt2_tt1;
                p = (2*tt2 + tt1 - 3*knotsites(i)) / tt2_tt1_sq;
                r = (2*knotsites(i) - tt2 - tt1) / (tt2_tt1_sq*tt2_tt1);
                for j = 1 : n
                    tmpx = X(j,knotdims(i));
                    if tmpx <= tt1
                        Xx(j) = 0;
                    elseif tmpx < tt2
                        % Xx(j) = p*(tmpx-tt1)^2 + r*(tmpx-tt1)^3; % explicit version
                        tmpx_tt1 = tmpx - tt1;
                        tmpx_tt1_sq = tmpx_tt1 * tmpx_tt1;
                        Xx(j) = p*tmpx_tt1_sq + r*tmpx_tt1_sq*tmpx_tt1;
                    else
                        Xx(j) = tmpx - knotsites(i);
                    end
                end
            else
                % p = (3*knotsites(i) - 2*tt1 - tt2) / (tt1 - tt2)^2; % explicit version
                % r = (tt1 + tt2 - 2*knotsites(i)) / (tt1 - tt2)^3; % explicit version
                tt1_tt2 = tt1 - tt2;
                tt1_tt2_sq = tt1_tt2 * tt1_tt2;
                p = (3*knotsites(i) - 2*tt1 - tt2) / tt1_tt2_sq;
                r = (tt1 + tt2 - 2*knotsites(i)) / (tt1_tt2_sq*tt1_tt2);
                for j = 1 : n
                    tmpx = X(j,knotdims(i));
                    if tmpx <= tt1
                        Xx(j) = knotsites(i) - tmpx;
                    elseif tmpx < tt2
                        % Xx(j) = p*(tmpx-tt2)^2 + r*(tmpx-tt2)^3; % explicit version
                        tmpx_tt2 = tmpx - tt2;
                        tmpx_tt2_sq = tmpx_tt2 * tmpx_tt2;
                        Xx(j) = p*tmpx_tt2_sq + r*tmpx_tt2_sq*tmpx_tt2;
                    else
                        Xx(j) = 0;
                    end
                end
            end
        end
        Xn = Xn .* Xx;
    end

end
return
