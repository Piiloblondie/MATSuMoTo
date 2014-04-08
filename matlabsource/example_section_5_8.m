function example_section_5_8
%EXAMPLE_SECTION_5_8.m executes the example of using MATSuMoTo for the 
%problems shown in Section 5.8 of the code manual
%--------------------------------------------------------------------------
%Copyright (c) 2014 by Juliane Mueller
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
data_file = 'datainput_hartman3'; %name of data file
maxeval = 300; %maximum number of allowed function evaluations
surogate_model = 'RBFtps'; %selected surrogate model type
sampling_technique = 'CANDglob'; %global randomized sampling strategy
initial_design = 'LHS'; %MATLAB's lhsdesign.m as initial design
number_startpoints = 15; %15 points in the initial experimental design
starting_point = []; %no user-specified points to be added to the initial design
NumberNewSamples = 1; %1 new point is selected in each iteration
[xbest, fbest] = MATSuMoTo(data_file,maxeval,surogate_model,sampling_technique,...
    initial_design,number_startpoints,starting_point,NumberNewSamples);

end %function