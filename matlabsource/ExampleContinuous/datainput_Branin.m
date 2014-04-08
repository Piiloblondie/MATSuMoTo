function Data = datainput_Branin
%DATAINPUT_BRANIN is an optimization test problem with only continuous 
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
Data.xlow=[-5,0]; %variable lower bounds
Data.xup=[10,15]; %variable upper bounds
Data.dim = 2; %problem dimension
Data.integer=[]; %indices of integer variables
Data.continuous=(1:2); %indices of continuous variables
%objective function
Data.objfunction= @(x)(x(:,2)-5.1*x(:,1).^2./(4*pi^2)+5*x(:,1)./pi-6).^2 + 10*(1-1/(8*pi))*cos(x(:,1))+10;
end %function