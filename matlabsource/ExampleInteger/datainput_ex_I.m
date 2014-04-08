function Data= datainput_ex_I
%DATAINPUT_EX_I is an optimization test problem with only integer variables
%adopted from http://www.aridolan.com/ga/gaa/MultiVarMin.html
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
Data.xlow=-100*ones(1,5);%variabe lower bounds
Data.xup=100*ones(1,5); %variabe upper bounds
Data.dim= 5; %problem dimension
Data.integer=(1:5); %indices of integer variables
Data.continuous = []; %indices of continuous variables
%define objective function
Data.objfunction=@(x) x(:,1).*sin(x(:,1)) + 1.7*x(:,2).*sin(x(:,1)) - ...
    1.5*x(:,3) - 0.1*x(:,4).*cos(x(:,4)+x(:,5)-x(:,1)) +...
    (0.2*x(:,5).^2-x(:,2)) - 1;
end %function