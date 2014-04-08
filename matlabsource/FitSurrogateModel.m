function [lambda,gamma,mmodel,beta, w_m] = FitSurrogateModel(Data)
%FITSURROGATEMODEL.m computes the parameters of the desired surrogate model
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
%Input: 
%Data - structure containing all information about the problem
%
%Output:
%lambda,gamma  - vectors, parameters of RBF model; lambda=gamma=[] if RBF model not used
%mmodel - structure, parameters for MARS model; mmodel=[] if MARS model not used
%beta - vector, parameters of polynomial regression model; beta=[] if polynomial model not used
%w_m - vector, contains the weights for the models in ensembles; w_m=[] if no ensemble model used
%--------------------------------------------------------------------------

%initialize parameters of surrogate models
lambda=[]; %for RBF
gamma=[]; %for RBF
mmodel=[]; %for MARS, http://www.cs.rtu.lv/jekabsons/Files/ARESLab.pdf
beta=[]; %for polynomial regression
w_m=[]; %weights for model ensembles

%---------RBF models----------
if strcmp(Data.surrogate_model,'RBFlin') %linear RBF model
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'linear'); 
elseif strcmp(Data.surrogate_model,'RBFtps') %thin plate spline RBF model
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'TPS'); 
elseif strcmp(Data.surrogate_model,'RBFcub')    %cubic RBF model
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'cubic'); 
%---------polynomial regression models----------    
elseif strcmp(Data.surrogate_model,'POLYlin')   %linear polynomial model 
    beta=POLY(Data.S,Data.Ymed,'lin');
elseif strcmp(Data.surrogate_model,'POLYquad')  %quadratic polynomial model        
    beta=POLY(Data.S,Data.Ymed,'quad'); 
elseif strcmp(Data.surrogate_model,'POLYquadr') %reduced quadratic polynomial model        
    beta=POLY(Data.S,Data.Ymed,'quadr');
elseif strcmp(Data.surrogate_model,'POLYcub')   %cubic polynomial model        
    beta=POLY(Data.S,Data.Ymed,'cub');
elseif strcmp(Data.surrogate_model,'POLYcubr')  %reduced cubic polynomial model               
    beta=POLY(Data.S,Data.Ymed,'cubr');
%---------MARS model----------    
elseif strcmp(Data.surrogate_model,'MARS') %multivariate adaptive regression spline (MARS)
    mmodel = aresbuild(Data.S,Data.Ymed); %uses ARESlab toolbox
    
%------------------------------
% add more models here if necessary. Alter accordingly other parts in the
% code (output, function handles,...)
%------------------------------

%---------Mixture models (2)----------    
elseif strcmp(Data.surrogate_model,'MIX_RcM') %ensemble of cubic RBF and MARS
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'cubic'); % cubic RBF
    mmodel = aresbuild(Data.S,Data.Ymed); %MARS
    w_m=DempsterFor2models(Data); %uses Dempster-Shafer theory to adjust model weights
elseif strcmp(Data.surrogate_model,'MIX_RcPc') %ensemble of cubic RBF and full cubic polynomial
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'cubic'); % cubic RBF
    beta=POLY(Data.S,Data.Ymed,'cub'); %cubic polynomial
    w_m=DempsterFor2models(Data); %uses Dempster-Shafer theory to adjust model weights
elseif strcmp(Data.surrogate_model,'MIX_RcPq') %ensemble of cubic RBF and quadratic polynomial 
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'cubic'); % cubic RBF
    beta=POLY(Data.S,Data.Ymed,'quad'); %quadratic polynomial
    w_m=DempsterFor2models(Data); %uses Dempster-Shafer theory to adjust model weights
elseif strcmp(Data.surrogate_model,'MIX_RcPqr') %ensemble of cubic RBF and reduced quadratic polynomial 
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'cubic'); % cubic RBF
    beta=POLY(Data.S,Data.Ymed,'quadr'); %reduced quadratic polynomial
    w_m=DempsterFor2models(Data); %uses Dempster-Shafer theory to adjust model weights
elseif strcmp(Data.surrogate_model,'MIX_RcPcr') %ensemble of cubic RBF and reduced cubic polynomial 
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'cubic'); % cubic RBF
    beta=POLY(Data.S,Data.Ymed,'cubr'); %reduced cubic polynomial
    w_m=DempsterFor2models(Data); %uses Dempster-Shafer theory to adjust model weights
    
% add more 2-model combinations here if necessary. Update function "DempsterFor2models" accordingly 

%---------Mixture models (3)----------        
elseif strcmp(Data.surrogate_model,'MIX_RcPcM') %ensemble of cubic RBF, cubic polynomial, and MARS
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'cubic'); %  cubic RBF
    beta=POLY(Data.S,Data.Ymed,'cub'); %cubic polynomial 
    mmodel = aresbuild(Data.S,Data.Ymed); %MARS
    w_m=DempsterFor3models(Data); %uses Dempster-Shafer theory to adjust model weights
    
% add more 3-model combinations here if necessary. Update function
% "DempsterFor3models" accordingly     
end

end%function