function Data= datainput_ex1221_MI
%DATAINPUT_EX1221_MI is an optimization test problem with continuous AND 
%integer variables; adopted from MINLPLib 
%http://www.gamsworld.org/minlp/minlplib
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
Data.xlow=zeros(1,5); %variable lower bounds
Data.xup=[10,10,10,1,1]; %variable upper bounds
Data.dim= 5; %problem dimension
Data.integer=(1:2); %indices of integer variables
Data.continuous=(3:5); %indices of continuous variables
Data.objfunction=@(x)-(- 2*x(:,1) - 3*x(:,2) - 1.5*x(:,3) - 2*x(:,4) + 0.5*x(:,5)); %objective function handle
end %function