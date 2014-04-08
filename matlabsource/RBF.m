function [lambda, gamma]=RBF(S,Y,flag)
%RBF.m computes the parameters of the radial basis function interpolant
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
%S - Sample site matrix (m x d), m=number of sample sites, d=dimension
%Y - Objective function values corresponding to points in S
%flag - string determining what type of RBF model should be used
%
%Output:
%lambda, gamma - vectors with RBF model parameters
%--------------------------------------------------------------------------

[m,n]=size(S); %determine how many points are in S and what is the dimension
P=[S,ones(m,1)]; %set up matrix augmented with vector of ones
R=pdist2(S,S); %compute pairwise dstances between points in S, pdist2 is MATLAB built-in function

if strcmp(flag,'cubic') %cubic RBF
    Phi= R.^3;
elseif strcmp(flag,'TPS') %thin plate spline RBF
    R(R==0)=1;
    Phi=R.^2.*log(R);
elseif strcmp(flag,'linear') %linear RBF
    Phi=R;
end

A=[Phi,P;P', zeros(n+1,n+1)]; %set up matrix for solving for parameters
RHS=[Y;zeros(n+1,1)]; %right hand side of linear system
params=A\RHS; %compute parameters
lambda=params(1:m); 
gamma=params(m+1:end); %parameters of polynomial tail

end %function