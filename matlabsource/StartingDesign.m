function S=StartingDesign(Data)
%STARTINGDESIGN.m generates initial experimental design 
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
%Data - structure, contains all problem information
%
%output:
%S - matrix with design points, (m x d), m=number of points, d = dimension
%--------------------------------------------------------------------------
%LHS - MATLAB's built-in Latin hypercube sampling (maximizes minimum
%distance)
%SLHD - symmetric Latin hypercube design
%CORNER - uses all or several corner points plus mid point 
%--------------------------------------------------------------------------

if strcmp(Data.initial_design,'LHS')%lhsdesign
    %create design, scaled to [0,1]
    S = lhsdesign(Data.number_startpoints,Data.dim,'criterion','maximin',...
            'iterations',20);
    %scale points to original hypercube size
    S=repmat(Data.xup-Data.xlow,Data.number_startpoints,1).*S+...
        repmat(Data.xlow,Data.number_startpoints,1);    
elseif strcmp(Data.initial_design,'SLHD')%symmetric Latin hypercube
    %create design, scaled to [0,1]
    S=SLHD(Data);
    %scale points to original hypercube size
    S=repmat(Data.xup-Data.xlow,size(S,1),1).*S...
        +repmat(Data.xlow,size(S,1),1);    
elseif strcmp(Data.initial_design,'CORNER') %corner point strategy
    S=cornerpoints(Data);
end
    
end %function