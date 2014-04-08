function xbest = MaximinD_c(Data,xselected)
%MAXIMIND_C.m finds the point that maximizes the minimum distance to 
%already evaluated and selected points -- for CONTINUOUS problems
%
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
%Input
%Data - structure, contains all information about the optimization problem
%xselected  - sample points already selected in this iteration
%
%Output
%xbest - points that maximizes the minimum distance to already evaluated
%and selected points
%--------------------------------------------------------------------------

fminconoptions=optimset('Display','off'); %do not show any algorithm output
myfun=@(x) find_shortest_dist(x,[Data.S;xselected]); %objective function
x0=Data.xlow+rand(1,Data.dim).*(Data.xup-Data.xlow); %initial  point for optimization
xbest = fmincon(myfun, x0, [],[],[],[],Data.xlow, Data.xup,[],...
            fminconoptions);
end %function

function d=find_shortest_dist(x,S) %objective function
%find shortest distance to already evaluated and selected points
[~,d] = knnsearch(S,x); %knnsearch is MATLAB built-in function
d=-d;
end %objective function