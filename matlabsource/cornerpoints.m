function samplesites=cornerpoints(Data)
%CORNERPOINTS.m generates an initial experimental design with points in 
%corners of the hyper rectangle defined by Data.xup and Data.xlow and adds 
%one point at the center of the hyper rectangle
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
%input:
%Data - structure with problem information
%
%output:
%samplesites - matrix with points in initial experimental design
%--------------------------------------------------------------------------

if Data.number_startpoints > 2^Data.dim+1
    error('There are not enough corners to create the requested number of starting points')
end

S=zeros(2^Data.dim,Data.dim); %initialize sample site matrix 
for ii = 1:Data.dim
    S(:,ii)=repmat([repmat(Data.xlow(ii),2^(Data.dim-ii),1);...
        repmat(Data.xup(ii),2^(Data.dim-ii),1)],2^(ii-1),1);
end


%if there are more corner points than desired number of starting  points,
%randomly select corner points
if size(S,1) >Data.number_startpoints-1 
    r=randperm(size(S,1));
    use=S(r(1:Data.number_startpoints-1),:); %select random corners
    samplesites=[use;Data.xlow+0.5*(Data.xup-Data.xlow)]; %add center point
else
    samplesites=[S;Data.xlow+0.5*(Data.xup-Data.xlow)];%add center point
end

end %function