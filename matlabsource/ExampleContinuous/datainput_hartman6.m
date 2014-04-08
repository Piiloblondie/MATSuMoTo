function Data=datainput_hartman6
%DATAINPUT_HARTMAN6 is an optimization test problem with only continuous 
%variables
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
%Input: None
%Output: Data - structure with optimization problem information
%--------------------------------------------------------------------------
%
Data.xlow =zeros(1,6);  % variable lower bounds
Data.xup=ones(1,6);     % variable upper bounds
Data.dim=6;             %problem dimension
Data.integer=[];        %indices of variables with integer constraints
Data.continuous=(1:6);  %indices of continuous variables 
Data.objfunction=@(x)myfun(x); %handle to objective function
end %function

function y=myfun(x) %objective function
c=[1,1.2,3,3.2]'; %c,A,b are data vectors
A=[10,3,17,3.5,1.7,8;0.05,10,17,0.1,8,14;3,3.5,1.7,10,17,8;17,8,0.05,10,0.1,14];
P=[0.1312,0.1696,0.5569,0.0124,0.8283,0.5886;...
    0.2329,0.4135,0.8307,0.3736,0.1004,0.9991;...
    0.2348,0.1451,0.3522,0.2883,0.3047,0.665;...
    0.4047,0.8828,0.8732,0.5743,0.1091,0.0381];
x=x(:)'; % make sure vector is row vector
y=-sum(c.*exp(-sum(A.*(repmat(x,4,1)-P).^2,2))); %compute objective function value
end %myfun