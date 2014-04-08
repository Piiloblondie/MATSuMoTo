function Data=datainput_nvs09alt_I
%DATAINPUT_NVS09ALT_I is an optimization test problem with only integer 
%variables adopted from http://www.gamsworld.org/minlp/minlplib/
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
Data.xlow=3*ones(1,10); %variable lower bounds
Data.xup=99*ones(1,10); %variable upper bounds
Data.dim=10; %problem dimension
Data.integer=(1:10); %indices of integer variables
Data.continuous=[]; %indices of continuous variables
%objective function
Data.objfunction=@(i) -(-((log(i(:,1) - 2)).^2 + (log(100 - i(:,1))).^2 + (log(i(:,2) - 2)).^2 +...
    (log(100 - i(:,2))).^2 + (log(i(:,3) - 2)).^2 + (log(100 - i(:,3))).^2 + (log(i(:,4) - 2)).^2 +...
    (log(100 - i(:,4))).^2 + (log(i(:,5) - 2)).^2 + (log(100 - i(:,5))).^2 + (log(i(:,6) - 2)).^2 +...
    (log(100 - i(:,6))).^2 + (log(i(:,7) - 2)).^2 + (log(100 - i(:,7))).^2 + (log(i(:,8) - 2)).^2 +...
    (log(100 - i(:,8))).^2 + (log(i(:,9) - 2)).^2 + (log(100 - i(:,9))).^2 + (log(i(:,10) - 2)).^2 +...
    (log(100 - i(:,10))).^2 - (i(:,1).*i(:,2).*i(:,3).*i(:,4).*i(:,5).*i(:,6).*i(:,7).*i(:,8).*i(:,9).*i(:,10))^0.2)) ; 
end %function