function trainParams = aresparams(maxFuncs, c, cubic, cubicFastLevel, ...
selfInteractions, maxInteractions, threshold, prune, useMinSpan, ...
useEndSpan, maxFinalFuncs)
% aresparams
% Creates a structure of ARES configuration parameter values for further
% use with aresbuild, arescv, or arescvc functions.
%
% Call:
%   trainParams = aresparams(maxFuncs, c, cubic, cubicFastLevel, ...
%                            selfInteractions, maxInteractions, ...
%                            threshold, prune, useMinSpan, useEndSpan, ...
%                            maxFinalFuncs)
%
% All the arguments of this function are optional. Empty values are also
% accepted (the corresponding default values will be used).
% 
% Input:
% For most applications, it can be expected that the most attention should
% be paid to the following parameters: maxFuncs, c, cubic, maxInteractions,
% and maybe threshold.
%   maxFuncs      : The maximal number of basis functions included in the
%                   model in the forward model building phase (before
%                   pruning in the backward phase). Includes the intercept
%                   term. (default value = 21). The recommended value for
%                   this parameter is two times the number of basis
%                   functions in the final model (Friedman 1991).
%   c             : Generalized Cross-Validation (GCV) penalty per knot.
%                   Theory suggests values in the range of about 2 to 4.
%                   Larger values for c will lead to fewer knots being
%                   placed (i.e., final models will be simpler). A value of
%                   0 penalizes only terms, not knots (can be useful e.g.,
%                   with lots of data and low noise). The recommended (and
%                   default) value for c is 3 (Friedman 1991). Note that if
%                   maxInteractions = 1 (additive modelling) then function
%                   aresbuild will recalculate c so that the actually used
%                   value is 2c/3 (recommended for additive modelling).
%   cubic         : Whether to use piecewise-cubic (true) or piecewise-
%                   linear (false) type of modelling (Friedman 1991). It is
%                   expected that the piecewise-cubic modelling will give
%                   higher predictive performance for smoother and less
%                   noisy data. (default value = true)
%   cubicFastLevel : In ARESLab, there are three types (levels) of
%                   piecewise-cubic modelling implemented. In level 0 cubic
%                   modelling for each candidate model is done in both
%                   phases of the technique (slow). In level 1 cubic
%                   modelling is done only in the backward phase (much
%                   faster). In level 2 cubic modelling is done after both
%                   phases only for the final model (fastest). The default
%                   and recommended level is 2. Levels 0 and 1 may bring
%                   extra precision in the modelling process however the
%                   results can actually also be worse. It is expected that
%                   the two much slower levels will mostly be not worth the
%                   waiting.
%   selfInteractions : The maximum degree of self interactions for any
%                   input variable. In ARESLab it can be larger than 1 only
%                   for piecewise-linear modelling. Usually the self
%                   interactions are never allowed. (default value = 1, no
%                   self interactions)
%   maxInteractions : The maximum degree of interactions between input
%                   variables. Set to 1 for additive modelling (i.e., no
%                   interactions). For maximal interactivity between the
%                   variables, set the parameter to d * selfInteractions,
%                   where d is the number of input variables – this way the
%                   modelling procedure will have the most freedom building
%                   a complex model. Typically only a low degree of
%                   interaction is allowed, but higher degrees can be used
%                   when the data warrants it. (default value = 1)
%   threshold     : One of the stopping criteria for the forward phase. The
%                   larger the value of threshold the potentially simpler
%                   models are generated (see reference manual for
%                   details). Default value = 1e-4. For noise-free data the
%                   value may be lowered.
%   prune         : Whether to perform the model pruning (the backward
%                   phase) or not. (default value = true)
%   useMinSpan    : Allows enabling, disabling, and manually tuning the
%                   protection agains runs of positive or negative error
%                   values between knots (see reference manual for
%                   details). (default and recommended value = -1 which
%                   corresponds to the automatic mode)
%   useEndSpan    : Allows enabling, disabling, and manually tuning the
%                   protection agains local variance of the estimates near
%                   the ends of data intervals (see reference manual for
%                   details). (default and recommended value = -1 which
%                   corresponds to the automatic mode)
%   maxFinalFuncs : Maximum number of basis functions (including the
%                   intercept term) in the pruned model. Use this (rather
%                   than the maxFuncs parameter) to enforce an upper bound
%                   on the final model size. (default value = Inf).
%
% Output:
%   trainParams   : A structure of training parameters for aresbuild
%                   function containing the provided values of the
%                   parameters (or default ones, if not provided).

% =========================================================================
% ARESLab: Adaptive Regression Splines toolbox for Matlab
% Author: Gints Jekabsons (gints.jekabsons@rtu.lv)
% URL: http://www.cs.rtu.lv/jekabsons/
%
% Copyright (C) 2009-2010  Gints Jekabsons
%
% =========================================================================

% Last update: December 11, 2009

if (nargin < 1) || (isempty(maxFuncs))
    trainParams.maxFuncs = 21;
else
    trainParams.maxFuncs = maxFuncs;
end

if (nargin < 2) || (isempty(c))
    trainParams.c = 3;
else
    trainParams.c = c;
end

if (nargin < 3) || (isempty(cubic))
    trainParams.cubic = true;
else
    trainParams.cubic = cubic;
end

if (nargin < 4) || (isempty(cubicFastLevel))
    trainParams.cubicFastLevel = 2;
else
    trainParams.cubicFastLevel = cubicFastLevel;
end

if (nargin < 5) || (isempty(selfInteractions))
    trainParams.selfInteractions = 1;
else
    trainParams.selfInteractions = selfInteractions;
end
if (trainParams.cubic) && (trainParams.selfInteractions > 1)
    trainParams.selfInteractions = 1;
    disp('Warning: trainParams.selfInteractions value reverted to 1 due to piecewise-cubic model setting.');
end

if (nargin < 6) || (isempty(maxInteractions))
    trainParams.maxInteractions = 1; % applicable maximum is d * trainParams.selfInteractions
else
    trainParams.maxInteractions = maxInteractions;
end

if (nargin < 7) || (isempty(threshold))
    trainParams.threshold = 1e-4;
else
    trainParams.threshold = threshold;
end

if (nargin < 8) || (isempty(prune))
    trainParams.prune = true;
else
    trainParams.prune = prune;
end

if (nargin < 9) || (isempty(useMinSpan))
    trainParams.useMinSpan = -1; % default = -1 = automatic
else
    if useMinSpan == 0
        trainParams.useMinSpan = 1; % 1 and 0 is the same here (no endspan)
    else
        trainParams.useMinSpan = useMinSpan;
    end
end

if (nargin < 10) || (isempty(useEndSpan))
    trainParams.useEndSpan = -1; % default = -1 = automatic
else
    if useEndSpan == 0
        trainParams.useEndSpan = 1; % 1 and 0 is the same here (no endspan)
    else
        trainParams.useEndSpan = useEndSpan;
    end
end

if (nargin < 11) || (isempty(maxFinalFuncs))
    trainParams.maxFinalFuncs = Inf;
else
    trainParams.maxFinalFuncs = maxFinalFuncs;
end

return
