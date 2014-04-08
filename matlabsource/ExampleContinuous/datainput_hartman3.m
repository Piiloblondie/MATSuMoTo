function Data=datainput_hartman3
%DATAINPUT_HARTMAN3 is an optimization test problem with only continuous 
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
Data.xlow =zeros(1,3); % variable lower bounds
Data.xup=ones(1,3);     % variable upper bounds
Data.dim = 3; %problem dimesnion
Data.integer = []; %indices of integer variables
Data.continuous = (1:3); %indices of continuous variables
Data.objfunction=@(x)myfun(x); %handle to objective function
end %function

function y=myfun(x) %objective function
c=[1,1.2,3,3.2]'; %c,A,b are data vectors
A=[3, 10, 30; 0.1,10,35; 3, 10, 30;0.1,10,35];
P=[0.3689    0.1170    0.2673
    0.4699    0.4387    0.7470
    0.1091    0.8732    0.5547
    0.0382    0.5743    0.8828];
x=x(:)'; % make sure vector is row vector
y=-sum(c.*exp(-sum(A.*(repmat(x,4,1)-P).^2,2))); %compute objective function value
end %myfun