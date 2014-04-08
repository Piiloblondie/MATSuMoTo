function Data = datainput_Shekel7
%DATAINPUT_SHEKEL7 is an optimization test problem with only continuous 
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
Data.xlow=[0,0,0,0]; %lower variable bounds
Data.xup=[10,10,10,10]; %upper variable bounds
Data.dim=4; %problem dimension
Data.integer=[]; %indices of integer variables
Data.continuous=(1:4); %indices of continuous variables
A=[4.0*ones(4,1), ones(4,1), 8*ones(4,1), 6*ones(4,1), [3,2,5,8,6,7;7,9,5,1,2,3.6; ...
    3,2,3,8,6,7;7,9,3,1,2,3.6]];
c=1/10*[1,2,2,4,4,6,3,7,5,5]; 
A=A(:,1:7);
c=c(1:7);
Data.objfunction=@(x)shekel(x,A,c); %handle to objective function
end %function

function y=shekel(x,A,c)
    x=x(:)';
    y=zeros(size(x,1),1);
    S1=zeros(size(A,2),1);
    for ii =1:size(x,1)
        for jj = 1:size(A,2)
            S1(jj,1)=1./(sum((x(ii,:)'-A(:,jj)).^2)+c(jj));
        end
        y(ii,1)=-sum(S1);
        S1=zeros(size(A,2),1);
    end
end %shekel