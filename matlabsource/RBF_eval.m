function Yest=RBF_eval(X,S,lambda,gamma,flag)
%RBF_EVAL.m uses radial basis function model to predict objective funtion 
%values
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
%input:
%X are points where function values should be calculated
%S are points where the function values are known
%lambda - parameter vector
%gamma contains the parameters of the polynomial tail
%flag is a string indicating which RBF to be used
%output: 
%the estimated function values at the points in X
%--------------------------------------------------------------------------

[mX,nX]=size(X); %dimensions of the points where function value should be predicted
[mS,nS]=size(S); %dimensions of sample sites taken so far 
if nX~=nS %check that both matrices are of the correct shape
    X=X';
    [mX,nX]=size(X);
end
R = pdist2(X,S); %compute pairwise dstances between points in X and S. pdist2 is MATLAB built-in function

if strcmp(flag,'cubic') %cubic RBF
    Phi= R.^3;
elseif strcmp(flag,'TPS') %thin plate spline RBF
    R(R==0)=1;
    Phi=R.^2.*log(R);
elseif strcmp(flag,'linear') %linear RBF
    Phi=R;
end
    
Yest1=Phi*lambda; %first part of response surface - weighted sum of distances
Yest2=[X,ones(mX,1)]*gamma; % polynomial tail of response surface
Yest=Yest1+Yest2; %predicted function value

end %functions