function InitialPoints = SLHD(Data)
% SLHD.m creates a symmetric Latin hypercube design. 
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
%Data - structure with all problem information 
%
%Output:
%InitialPoints - points in the starting design
%--------------------------------------------------------------------------

delta = (1/Data.number_startpoints)*ones(1,Data.dim);

X = zeros(Data.number_startpoints,Data.dim);
for j = 1:Data.dim
    for i = 1:Data.number_startpoints
        X(i,j) = ((2*i-1)/2)*delta(j);
    end
end

P = zeros(Data.number_startpoints,Data.dim);
P(:,1) = (1:Data.number_startpoints)';
if (mod(Data.number_startpoints,2) == 0)
   k = Data.number_startpoints/2;
else
   k = (Data.number_startpoints-1)/2;
   P(k+1,:) = (k+1)*ones(1,Data.dim);
end

for j = 2:Data.dim
   P(1:k,j) = randperm(k)';
   for i = 1:k
      if (rand(1) <= 0.5)
         P(Data.number_startpoints+1-i,j) = Data.number_startpoints+1-P(i,j);
      else
         P(Data.number_startpoints+1-i,j) = P(i,j);
         P(i,j) = Data.number_startpoints+1-P(i,j);
      end
   end
end

InitialPoints = zeros(Data.number_startpoints,Data.dim);
for j = 1:Data.dim
    for i = 1:Data.number_startpoints
        InitialPoints(i,j) = X(P(i,j),j);
    end
end

end%function