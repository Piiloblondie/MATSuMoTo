function xselected = SurfMin(Data,lambda,gamma, mmodel,beta, w_m)
%SURFMIN.m uses minimum point of response surface to select the next sample
%point in each iteration    
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
%Input
%Data - structure, contains all information about the optimization problem
%lambda,gamma  - vectors, parameters of RBF model; lambda=gamma=[] if RBF model not used
%mmodel - structure, parameters for MARS model; mmodel=[] if MARS model not used
%beta - vector, parameters of polynomial regression model; beta=[] if polynomial model not used
%w_m - vector, contains the weights for the models in mixtures; w_m=[] if no ensemble model used
%
%Output
%xselected - points for doing the next expensive objective function
%evaluations
%--------------------------------------------------------------------------


%determine which respose surface to minimize
%------------RBF models ----------------
if strcmp(Data.surrogate_model,'RBFlin') %linear RBF model
    myfun=@(x)RBF_eval(x,Data.S,lambda,gamma,'linear'); 
elseif strcmp(Data.surrogate_model,'RBFtps') %thin-plate spline RBF model
    myfun=@(x)RBF_eval(x,Data.S,lambda,gamma,'TPS'); 
elseif strcmp(Data.surrogate_model,'RBFcub')  % cubic RBF model  
    myfun=@(x)RBF_eval(x,Data.S,lambda,gamma,'cubic'); 
%------------Polynomial models --------------     
elseif strcmp(Data.surrogate_model,'POLYlin')  %linear polynomial model
    myfun=@(x)POLY_eval(x,beta,'lin');
elseif strcmp(Data.surrogate_model,'POLYquad')  %quadratic polynomial model       
    myfun=@(x)POLY_eval(x,beta,'quad'); 
elseif strcmp(Data.surrogate_model,'POLYquadr') %reduced quadratic polynomial model
    myfun=@(x)POLY_eval(x,beta,'quadr');
elseif strcmp(Data.surrogate_model,'POLYcub')  %cubic polynomial model       
    myfun=@(x)POLY_eval(x,beta,'cub'); 
elseif strcmp(Data.surrogate_model,'POLYcubr')  %reduced cubic polynomial model      
    myfun=@(x)POLY_eval(x,beta,'cubr');
%------------MARS models --------------     
elseif strcmp(Data.surrogate_model,'MARS') %multivariate adaptive regression splines
    myfun=@(x)arespredict(mmodel,x);
%------------Ensemble models (2)--------------     
elseif strcmp(Data.surrogate_model,'MIX_RcM') %Ensemble: cubic RBF and MARS
    myfun1=@(x)RBF_eval(x,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model
    myfun2=@(x)arespredict(mmodel,x); %objective function value prediction with MARS model
    myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x); %weighted objective function value prediction   
elseif strcmp(Data.surrogate_model,'MIX_RcPc')  %Ensemble: cubic RBF and full cubic polynomial
    myfun1=@(x)RBF_eval(x,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model 
    myfun2=@(x)POLY_eval(x,beta,'cub');%objective function value prediction with full cubic polynomial 
    myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x); %weighted objective function value prediction            
elseif strcmp(Data.surrogate_model,'MIX_RcPq') %Ensemble: cubic RBF and full quadratic polynomial
    myfun1=@(x)RBF_eval(x,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model 
    myfun2=@(x)POLY_eval(x,beta,'quad');%objective function value prediction with quadratic polynomial
    myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x); %weighted objective function value prediction    
elseif strcmp(Data.surrogate_model,'MIX_RcPcr') %Ensemble: cubic RBF and reduced cubic polynomial
    myfun1=@(x)RBF_eval(x,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model 
    myfun2=@(x)POLY_eval(x,beta,'cubr');   %objective function value prediction with reduced cubic polynomial
    myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x); %weighted objective function value prediction  
elseif strcmp(Data.surrogate_model,'MIX_RcPqr') %Ensemble: cubic RBF and reduced quadratic polynomial
    myfun1=@(x)RBF_eval(x,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model 
    myfun2=@(x)POLY_eval(x,beta,'quadr');   %objective function value prediction with reduced quadratic polynomial
    myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x); %weighted objective function value prediction  
%add other 2-model mixtures here

%------------Mixture models (3)--------------         
elseif strcmp(Data.surrogate_model,'MIX_RcPcM') %Ensemble: cubic RBF, full cubic polynomial, and MARS       
    myfun1=@(x)RBF_eval(x,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model
    myfun2=@(x)POLY_eval(x,beta,'cub');%objective function value prediction with full cubic polynomial 
    myfun3=@(x)arespredict(mmodel,x); %objective function value prediction with MARS model
    myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x)+ w_m(3)*myfun3(x); %weighted objective function value prediction  
end

if isempty(Data.integer) %continuous problem, use fmincon
    %fmincon options for finding minimum of response surface
    fminconoptions=optimset('Display','off');
    xselected=[]; %initialize matrix with new sample points
    for ii = 1:Data.NumberNewSamples 
        x0=Data.xlow+rand(1,Data.dim).*(Data.xup-Data.xlow); %initial  point for optimization
        xbest=fmincon(myfun, x0, [],[],[],[],Data.xlow, Data.xup,[],...
            fminconoptions);
        seltry=0;
        %in case selected point is already contained in S or already
        %selected: try to find a different point and if without success,
        %use the point that maximizes the minimum distance to the already
        %sampled points
        while any(sqrt(sum((repmat(xbest,size([Data.S;xselected],1),1)-[Data.S;xselected]).^2,2))<Data.tolerance)
            seltry=seltry+1;
            x0=Data.xlow+rand(1,Data.dim).*(Data.xup-Data.xlow); %initial point for optimization
            xbest=fmincon(myfun, x0, [],[],[],[],Data.xlow, Data.xup,[],...
                fminconoptions);
            if seltry>3 %3 trials couldn't find a new point, select the point that maximizes the minimal distance to all other points
                xbest = MaximinD_c(Data,xselected);
                while any(sqrt(sum((repmat(xbest,size([Data.S;xselected],1),1)-[Data.S;xselected]).^2,2))<Data.tolerance)
                    xbest = Data.xlow+(Data.xup-Data.xlow).*rand(1,Data.dim);
                end
            end
        end
        xselected=[xselected;xbest];
    end
else
    %GA for mixed integer
    GAoptions = gaoptimset('Display','off');
    xselected=[];
    for ii = 1:Data.NumberNewSamples
        xbest = ga(myfun,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer,GAoptions);
        seltry=0;
        %in case selected point is already contained in S or already
        %selected: try to find a different point and if without success,
        %use the point that maximizes the minimum distance to the already
        %sampled points
        while any(sqrt(sum((repmat(xbest,size([Data.S;xselected],1),1)-[Data.S;xselected]).^2,2))<Data.tolerance)
            seltry=seltry+1;
            xbest = ga(myfun,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer,GAoptions);
            if seltry>3 %3 trials couldnt find a new point, select the point that maximizes the minimal distance to all other points
                xbest=MaximinD_i(Data,xselected);
                while any(sqrt(sum((repmat(xbest,size([Data.S;xselected],1),1)-[Data.S;xselected]).^2,2))<Data.tolerance)
                    xbest = Data.xlow+(Data.xup-Data.xlow).*rand(1,Data.dim);
                    xbest(Data.integer)=round(xbest(Data.integer));
                end
            end
        end
        xselected=[xselected;xbest];
    end
end

end%function