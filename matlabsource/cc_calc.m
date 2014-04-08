function corcoef=cc_calc(Yo,Yp)
%CORCOEF.m calculates correlation coefficiants of model response and true 
%function value of cross validation
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
%Input: 
%Yo - observed function values
%Yp - predicted function values
%Output:
%corcoef - correlation coefficients
%--------------------------------------------------------------------------

corcoef=zeros(size(Yp,2),1);
for ii = 1:size(Yp,2)
    CCmatrix=corrcoef(Yo, Yp(:,ii));
    corcoef(ii)=CCmatrix(1,2);
end
end%function