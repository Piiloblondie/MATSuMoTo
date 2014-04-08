function CandPoint = Perturbation(Data, NCandPoint, sigma_stdev,P)
%PERTURBATION.m generates candidate points for the next sample site by 
%perturbing vaiables of best point found so far. This routine is for
%CONTINUOUS problems
%If 'CANDglob' sampling, generate points also by uniformly selecting points
%from whole variable domain
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
%Input:
% Data - structure with all problem data
% NCandPoint - number of candidate points generated for each group
% sigma_stdev - perturbation rate for continuous variables
% P - probability with which each variable is perturbed
%Output:
% CandPoint - candidate points
%--------------------------------------------------------------------------

%Group 1: perturbations of best point found so far
CandPoint_G1=repmat(Data.xbest, NCandPoint,1);
for ii = 1:NCandPoint %generate NCandPoint candidate points
    rp=rand(1,Data.dim);
    pertID=zeros(1,Data.dim);
    pertID(rp<P)=1; %find the variables of xbest to perturb
    if sum(pertID)==0 %no variables are to be perturbed -> select one variable at random
        pID=randperm(Data.dim);
        pertID(pID(1))=1;
    end
    for jj = 1:Data.dim %go through each dimension
        if pertID(jj)==1 %check if perturbation should be done
            p=randperm(length(sigma_stdev)); %randomly pick the rate for pertubation
            CandPoint_G1(ii,jj) = max(Data.xlow(jj), ...
                min(CandPoint_G1(ii,jj)+sigma_stdev(p(1))*randn(1),...
                Data.xup(jj)));
        end
    end
end
    
if strcmp(Data.sampling_technique,'CANDglob')
    %Group 2: uniformly select variabe values from whole variable domain
    CandPoint_G2=repmat(Data.xlow,NCandPoint,1) + ...
        rand(NCandPoint,Data.dim).* repmat(Data.xup-Data.xlow,NCandPoint,1);
    CandPoint=[CandPoint_G1;CandPoint_G2];
    clear CandPoint_G1 CandPoint_G2;
else
    CandPoint=CandPoint_G1;
    clear CandPoint_G1
end

end %function