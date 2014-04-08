function TestDriver
%TESTDRIVER.m is a test driver to run MATSuMoTo for a test problem without
%user input.
%Click on RUN to execute the program.
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
%
disp('Testing continuous problem')
MATSuMoTo('datainput_hartman3',150,'RBFcub','CANDglob','SLHD',20,[],1);
disp('Continuous problem finished successfully')

disp('Testing integer problem')
MATSuMoTo('datainput_convex_I',100,'RBFcub','CANDglob','SLHD',20,[],1);
disp('Integer problem finished successfully')

disp('Testing mixed-integer problem')
MATSuMoTo('datainput_ex_MI',200,'RBFcub','CANDglob','SLHD',20,[],1);
disp('Mixed-integer problem finished successfully')