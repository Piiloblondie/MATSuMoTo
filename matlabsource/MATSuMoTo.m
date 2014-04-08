function [xbest,fbest]=MATSuMoTo(data_file,varargin)
%MATSUMOTO - MATLAB Surrogate Model Toolbox for unconstrained continuous, 
%integer, and constrained mixed-integer global optimization problems. 
%MATSUMOTO attempts to solve problems of the form
%          min F(X)
%           X    
%          s.t. LB <= X <= UB
%
%--------------------------------------------------------------------------
%Copyright (c) 2013 by Juliane Mueller
%
%This file is part of MATSuMoTo.m - the MATLAB Surrogate Model Toolbox
%MATSuMoTo is free software: you can redistribute it and/or modify it under
%the terms of the GNU General Public License as published by the Free 
%Software Foundation, either version 3 of the License, or (at your option) 
%any later version.
%
%MATSuMoTo is distributed in the hope that it will be useful, but WITHOUT 
%ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
%FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
%more details.
%
%You should have received a copy of the GNU General Public License along 
%with MATSuMoTo.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
% Function call:
% MATSuMoTo(data_file,maxeval,surogate_model,sampling_technique,...
%    initial_design, number_startpoints,starting_point,NumberNewSamples);
% where 
% data_file - mandatory, string with name of data file containing 
%               optimization problem definition (see code documentation)
% maxeval - optional, positive integer defining maximum number of allowed 
%               function evaluations (default 20*dimension)
% surrogate_model- optional, string defining which surrogate model to  use 
%               (default 'RBFcub')
% sampling_technique - optional, string defining the technique for 
%               selecting the next sample site (default 'CANDglob')
% initial_design - optional, string defining the initial experimental 
%               design technique (default 'LHS')
% number_startpoints - optional, positive integer defining the number of 
%               initial starting points (default 2*(dimension+1))
% starting_point - optional, matrix with specific points to be included in 
%               initial experimental design (default [])
% NumberNewSamples - optional, positive integer determining the number of 
%               points to be evaluated in every iteration (default 1)
%
%The user can choose between different toolbox options for
%   (A) the surrogate model: 'RBFcub', 'RBFtps', 'RBFlin', 'POLYlin', 
%               'POLYquad', 'POLYquadr', 'POLYcub', 'POLYcubr', 'MARS',
%               'MIX_RcM', 'MIX_RcPc', 'MIX_RcPq', 'MIX_RcPcM', 
%               'MIX_RcPqr', 'MIX_RcPcr' 
%   RBFcub - cubic RBF (default)
%   RBFtps - thin-plate spline RBF
%   RBFlin - linear RBF
%   POLYlin - linear regression polynomial
%   POLYquad - quadratic regression polynomial
%   POLYquadr - reduced quadratic regression polynomial
%   POLYcub - cubic regression polynomial
%   POLYcubr - reduced cubic regression polynomial
%   MARS - multivariate adaptive regression spline (uses ARESlab toolbox)
%   MIX_RcM - ensemble of cubic RBF and MARS
%   MIX_RcPc - ensemble of cubic RBF and full cubic regression polynomial
%   MIX_RcPq - ensemble of cubic RBF and full quadratic regression polynomial
%   MIX_RcPcM - ensemble of cubic RBF, full cubic regression polynomial, and MARS
%   MIX_RcPqr - ensemble of cubic RBF and reduced quadratic regression polynomial
%   MIX_RcPcr - ensemble of cubic RBF and reduced cubic regression polynomial
%   (B) the sampling strategy: 'CANDglob', 'CANDloc', 'SurfMin'
%   CANDglob - global randomized sampling
%   CANDloc - local randomized sampling
%   SurfMin - minimum point of surrogate model sampling
%   (C) the initial experimental design: 'LHS', 'SLHD', 'CORNER'
%   LHS - MATLAB's built-in Latin hypercube sampling
%   SLHD - symmetric Latin hypercube sampling
%   CORNER - uses (subset of) corner points and midpoint of hypercube
% 
% Example:
% MATSuMoTo('datainput_hartman3', 300, 'RBFtps', 'CANDglob', 'LHS', 15, [], 1);
% This example attempts to find the minimum of the three-dimensional 
%   Hartmann function. The maximum number of function evaluations is set to
%   300, the surrogate model to be used is the thin-plate spline radial 
%   basis function. The global randomized sampling strategy is used, and 
%   the initial experimental design is built with a Latin hypercube design
%   with 15 points. Starting points are not given ([]) and only one point 
%   is selected in every iteration.
%
% The results are saved in the file Results.mat in the current working
% directory. To access, type 
% load Results.mat
% into the MATLAB command window.
%--------------------------------------------------------------------------
%
% MATSuMoTo is based on various published papers that should be cited
% when using the toolbox:
% 1) J. Mueller and R. Piche, 2011, "Mixture Surrogate Models Based on 
%    Dempster-Shafer Theory for Global Optimization Problems", Journal of 
%    Global Optimization, 51:79-104
% 2) J. Mueller, C.A. Shoemaker, and R. Piche, 2013, "SO-MI: A Surrogate 
%    Model Algorithm for Computationally Expensive Nonlinear Mixed-Integer 
%    Black-Box Global Optimization Problems", Computers and Operations 
%    Research, 40(5):1383-1400
% 3) J. Mueller, C.A. Shoemaker, and  R. Piche, 2013, "SO-I: A Surrogate 
%    Model Algorithm for Expensive Nonlinear Integer Programming Problems 
%    Including Global Optimization Applications", Journal of Global 
%    Optimization, DOI 10.1007/s10898-013-0101-y
% 4) J. Mueller and C.A. Shoemaker, 2014, "Influence of Ensemble Surrogate 
%    Models and Sampling Strategy on the Solution Quality of Algorithms for
%    Computationally Expensive Black-Box Global Optimization Problems",  
%    Journal of Global Optimization, DOI 10.1007/s10898-014-0184-0
%--------------------------------------------------------------------------
%
timeStart=tic; %record time of algorithm

%% ---------- START INPUT CHECK -------------------
Data = InputCheck(data_file,varargin{:});
Data.Problem=data_file; %save the name of the problem to the structure
%% ---------- END INPUT CHECK ---------------------

%% ------ START OPIMIZATION -----
%optimization phase - distinguish between continuous/integer/mixed-integer
if isempty(Data.continuous) && ~isempty(Data.integer) %purely integer problem
    Solution= Optimization_integer(Data);    
elseif ~isempty(Data.continuous) && isempty(Data.integer) %purely continuous problem 
    Solution= Optimization_continuous(Data);   
elseif ~isempty(Data.continuous) && ~isempty(Data.integer) %mixed-integer problem   
    Solution= Optimization_mixedinteger(Data);    
end
%% ---update elements of Data strructure with results from the optimization
Data.S = Solution.S; %sample points
Data.Y = Solution.Y; %function values
Data.fbest = Solution.fbest; %best function value found in each optimization trial
Data.xbest = Solution.xbest; %best point found in each optimization trial
Data.fevaltime = Solution.fevaltime; %time for each function evaluation
Data.EvalsEachTrial = Solution.EvalsEachTrial; % number of function evaluations in each trial
Data.trial = Solution.trial; %number of trials in each optimization run
Data.TotalTime=toc(timeStart); %total algorithm run time

%% ------ SAVE RESULTS --------
save('Results','Data'); 

%% ---------Progress plot --------
%need to sort the objective function values in decreasing fashion  
plotY=zeros(length(Data.Y),1);
plotY(1)=Data.Y(1);
for ii = 2:length(Data.Y)
    if Data.Y(ii) <plotY(ii-1)
        plotY(ii)=Data.Y(ii);
    else
        plotY(ii) = plotY(ii-1);
    end
end
figure
plot((1:length(plotY)),plotY) %might want to use semilogy(...)
xlabel('Number of function evaluations')
ylabel('Objective function value')
title('Progress plot for objective function value')
%some statistics to print out at the end of the optimization
[bestFval,IDbest]=min(Data.fbest); %overall best function value
bestPoint = Data.xbest(IDbest,:); %overall best sample point
FF=repmat('%f ',1,Data.dim); 
fprintf('My best feasible function value is: %f \n',bestFval)
fprintf(['My best feasible point is: [',num2str(FF),'] \n'],bestPoint)
fprintf('Your output is saved in the file Results.mat in the current working directory. \n Thank you. Come again.\n')
xbest=bestPoint;
fbest=bestFval;
if matlabpool('size') > 0 %close the matlabpool if necessary
    matlabpool close
end
end %function MATSuMoTo