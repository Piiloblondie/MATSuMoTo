function Data= datainput_Rastrigin12
%DATAINPUT_RASTRIGIN12 is an optimization test problem with only continuous 
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
n=12;
Data.xlow=-ones(1,n); %lower variable bounds
Data.xup=3*ones(1,n); %upper variable bounds
Data.dim = n; %problem dimension
Data.integer = []; %integer variable indices
Data.continuous = (1:n); %continuous variable indices
Data.objfunction=@(x)myfun(x);  %objective function handle
end %function

function y=myfun(x)
x=x(:)';
y=sum((x.^2) - cos(2*pi*x),2);
end %myfun