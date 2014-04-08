function Data = datainput_Ackley15
%DATAINPUT_ACKLEY15 is an optimization test problem with only continuous 
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
Data.xlow=-15*ones(1,15); %lower variable bounds
Data.xup=30*ones(1,15); %variable upper bounds
Data.dim=15; %problem dimension
Data.integer =[]; %indices of integer variables
Data.continuous = (1:15); %indices of continuous variables
%objective function
Data.objfunction=@(x)-20*exp(-0.2*sqrt(sum(x.^2,2)/Data.dim)) - exp(sum(cos(2*pi*x),2)/Data.dim);
end %function