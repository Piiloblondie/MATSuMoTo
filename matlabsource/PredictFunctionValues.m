function CandValue=PredictFunctionValues(Data,CandPoint,lambda,gamma,...
    beta,mmodel,w_m)
%PREDICTFUNCTIONVALUES.m predicts the objective function values using the 
%desired surrogate model
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
%input
%Data - structure containing all information about the optimization problem
%CandPoint - matrix with candidates for next sample point. For these points
%the objective function values will be predicted
%lambda,gamma - parameters of RBF model, lambda=gamma=[] if RBF model not
%used
%beta - parameters of polyniomial regression model, beta = [] if
%polynomial regression model not used
%mmodel - structure with parameters of MARS model, mmodel = [] if MARS
%model not used (uses ARESlab toolbox)
%w_m - weights for models in mixtures, w_m = [] if ensemble model not used
%
%Output
%CandValue - predicted objective function values for points in CandPoint
%--------------------------------------------------------------------------


%----------RBF models---------
if strcmp(Data.surrogate_model,'RBFlin') %linear RBF model
    CandValue=RBF_eval(CandPoint,Data.S,lambda,gamma,'linear'); 
elseif strcmp(Data.surrogate_model,'RBFtps') %thin-plate spline RBF model
    CandValue=RBF_eval(CandPoint,Data.S,lambda,gamma,'TPS'); 
elseif strcmp(Data.surrogate_model,'RBFcub')  % cubic RBF model  
    CandValue=RBF_eval(CandPoint,Data.S,lambda,gamma,'cubic'); 
%----------Polynomial models---------
elseif strcmp(Data.surrogate_model,'POLYlin')  %linear polynomial
    CandValue=POLY_eval(CandPoint,beta,'lin');
elseif strcmp(Data.surrogate_model,'POLYquad') %full quadratic polynomial        
    CandValue=POLY_eval(CandPoint,beta,'quad'); 
elseif strcmp(Data.surrogate_model,'POLYquadr') %reduced quadratic polynomial        
    CandValue=POLY_eval(CandPoint,beta,'quadr');
elseif strcmp(Data.surrogate_model,'POLYcub')  %full cubic polynomial               
    CandValue=POLY_eval(CandPoint,beta,'cub'); 
elseif strcmp(Data.surrogate_model,'POLYcubr')  %reduced cubic polynomial                     
    CandValue=POLY_eval(CandPoint,beta,'cubr');
%----------MARS model---------    
elseif strcmp(Data.surrogate_model,'MARS') %Multivariate adaptive regression spline
    CandValue=arespredict(mmodel,CandPoint); %uses ARESlab toolbox
%----------Ensemble models (2)---------    
elseif strcmp(Data.surrogate_model,'MIX_RcM')  %Ensemble: cubic RBF & MARS
    CandValue_Rc=RBF_eval(CandPoint,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model 
    CandValue_M=arespredict(mmodel,CandPoint); % objective function value prediction with MARS model 
    CandValue= w_m(1)*CandValue_Rc+w_m(2)*CandValue_M; %weighted objective function value prediction
elseif strcmp(Data.surrogate_model,'MIX_RcPc')  %Ensemble: cubic RBF & full cubic polynomial
    CandValue_Rc=RBF_eval(CandPoint,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model 
    CandValue_Pc=POLY_eval(CandPoint,beta,'cub');%objective function value prediction with full cubic polynomial 
    CandValue= w_m(1)*CandValue_Rc+w_m(2)*CandValue_Pc; %weighted objective function value prediction            
elseif strcmp(Data.surrogate_model,'MIX_RcPq')  %Ensemble: cubic RBF & quadratic polynomial
    CandValue_Rc=RBF_eval(CandPoint,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model 
    CandValue_Pq=POLY_eval(CandPoint,beta,'quad');%objective function value prediction with quadratic polynomial
    CandValue= w_m(1)*CandValue_Rc+w_m(2)*CandValue_Pq; %weighted objective function value prediction            
elseif strcmp(Data.surrogate_model,'MIX_RcPcr')  %Ensemble: cubic RBF  & reduced cubic polynomial
    CandValue_RBF = RBF_eval(CandPoint,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model
    CandValue_Pcr=POLY_eval(CandPoint,beta,'cubr'); %objective function value prediction with reduced cubic polynomial 
    CandValue= w_m(1)*CandValue_RBF+w_m(2)*CandValue_Pcr;   %weighted objective function value prediction
elseif strcmp(Data.surrogate_model,'MIX_RcPqr')  %Ensemble: cubic RBF & reduced quadratic polynomial
    CandValue_RBF = RBF_eval(CandPoint,Data.S,lambda,gamma,'cubic');  % objective function value prediction with cubic RBF model
    CandValue_Pqr=POLY_eval(CandPoint,beta,'quadr'); %objective function value prediction with reduced quadratic polynomial
    CandValue= w_m(1)*CandValue_RBF+w_m(2)*CandValue_Pqr;  %weighted objective function value prediction

    
%add more surrogate models here if desired
%----------Ensemble models (3)---------        
elseif strcmp(Data.surrogate_model,'MIX_RcPcM')   %Ensemble: cubic RBF & full cubic polynomial & MARS        
    CandValue_Rc=RBF_eval(CandPoint,Data.S,lambda,gamma,'cubic');  % objective function value prediction with cubic RBF model 
    CandValue_Pcr=POLY_eval(CandPoint,beta,'cub');%objective function value prediction with full cubic polynomial 
    CandValue_M=arespredict(mmodel,CandPoint);  % objective function value prediction with MARS model 
    CandValue= w_m(1)*CandValue_Rc+w_m(2)*CandValue_Pcr + w_m(3)*CandValue_M; %weighted objective function value prediction
end

end%function