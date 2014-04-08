function [t1 t2 dif] = findsideknots(model, knotdim, knotsite, d, minX, maxX, old_t1, old_t2)
% Recomputes side knot placements for all the basis functions in the model
% (in either the knotdim dimension or all the dimensions) while taking into
% account that a new basis function is to be added and its side knots must
% not disrupt the integrity of the whole side knot placement. The function
% is used only for piecewise-cubic modelling.
% The function assumes that there is no self-interactions in any basis
% function.
% t1 is a matrix of knot sites for the knots on the left of the central
% knot. t2 is a matrix of knot sites for the knots on the right of the
% central knot. dif is a boolean vector with 1's for the basis functions
% for which at least one knot has moved, and 0's for other basis functions.

% =========================================================================
% ARESLab: Adaptive Regression Splines toolbox for Matlab
% Author: Gints Jekabsons (gints.jekabsons@rtu.lv)
% URL: http://www.cs.rtu.lv/jekabsons/
%
% Copyright (C) 2009-2010  Gints Jekabsons
%
% =========================================================================

% Last update: November 9, 2009

if isempty(model.knotdims)
    if isempty(knotdim)
        t1 = [];
        t2 = [];
        dif = [];
    else
        t1 = inf(1,d);
        t2 = inf(1,d);
        for di = 1 : d
            temp = find(knotdim == di);
            if ~isempty(temp)
                t1(di) = (minX(di) + knotsite(temp)) / 2;
                t2(di) = (knotsite(temp) + maxX(di)) / 2;
            end
        end
        dif = true;
    end
else
    if ~isempty(knotdim)
        t1 = old_t1;
        t1(end+1,:) = Inf;
        t2 = old_t2;
        t2(end+1,:) = Inf;
        dims = knotdim; % take dimension numbers from the new knot
    else
        t1 = inf(length(model.knotdims),d);
        t2 = inf(length(model.knotdims),d);
        dims = 1 : d;
    end
    % "ind" is indexes of functions which are in the di dimension
    % "knot" is knotsites of those functions in the di dimension
    for di = dims
        % find all used knotsites in the di dimension
        ind = zeros(length(model.knotdims)+1,1);
        knot = ind;
        count = 0;
        for i = 1 : length(model.knotdims)
            temp = model.knotdims{i} == di;
            if any(temp)
                count = count + 1;
                ind(count) = i;
                knot(count) = model.knotsites{i}(temp);
            end
        end
        if ~isempty(knotdim)
            temp = knotdim == di;
            if any(temp)
                count = count + 1;
                % here ind(count) is already 0
                knot(count) = knotsite(temp);
            end
        end
        ind = ind(1:count);
        knot = knot(1:count);
        % sort the arrays by knot place
        [knot ind2] = sort(knot);
        ind = ind(ind2);
        % calculate and store t1, t2
        % (it is possible that two or more functions have the same central
        % knot and so the same side knots)
        i = 1;
        while (i <= length(knot))
            if i == 1
                temp_t1 = (minX(di) + knot(i)) / 2;
            else
                temp_t1 = (knot(i-1) + knot(i)) / 2;
            end
            temp = find(knot == knot(i));
            for j = 1 : length(temp)
                if ind(temp(j)) ~= 0
                    t1(ind(temp(j)),di) = temp_t1;
                else
                    t1(end,di) = temp_t1;
                end
            end
            if knot(temp(end)) == knot(end)
                temp_t2 = (knot(i) + maxX(di)) / 2;
            else
                temp_t2 = (knot(i) + knot(temp(end)+1)) / 2;
            end
            for j = 1 : length(temp)
                if ind(temp(j)) ~= 0
                    t2(ind(temp(j)),di) = temp_t2;
                else
                    t2(end,di) = temp_t2;
                end
            end
            i = i + length(temp);
        end
    end
    % find the difference between old t1,t2 and new t1,t2
    if ~isempty(old_t1)
        if isempty(knotdim)
            dif = any((t1 ~= old_t1) | (t2 ~= old_t2), 2);
        else
            dif = any((t1(1:end-1,:) ~= old_t1) | (t2(1:end-1,:) ~= old_t2), 2);
        end
    else
        dif = [];
    end
end
return
