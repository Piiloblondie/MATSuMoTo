function Data= datainput_convex_I
%DATAINPUT_CONVEX_I is an optimization test problem with only integer 
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
Data.xlow=-10*ones(1,8); %variable lower bounds
Data.xup=10*ones(1,8); %variable upper bounds
Data.dim=8; %problem dimension
Data.integer=(1:8); %indices of integer variables
Data.continuous = []; %indices of continuous variables
%define objective function
Data.objfunction=@(x) 3.1*x(:,1).^2 + 7.6* x(:,2).^2 +6.9*x(:,3).^2 +0.004*x(:,4).^2 +...
    +19*x(:,5).^2 +3*x(:,6).^2 +x(:,7).^2  +4*x(:,8).^2 ;
end %function