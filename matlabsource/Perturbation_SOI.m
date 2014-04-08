function CandPoint = Perturbation_SOI(Data,NCandidates,stdev_int,P)
%PERTURBATION_SOI.m generates candidate points for the next sample site by 
%perturbing vaiables of best point found so far. This routine is for
%PURE INTEGER problems
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
% Data - structure, contains information about optimization problem
% NCandidates - integer, number of candidate points contained in each group
% stdev_int - integer vector, perturbation ranges
% P - perturbation probability
%
%Output:
%CandPoint - matrix with candidate points
%--------------------------------------------------------------------------

%set candidate points same as xbest, and then perturb only integer
%variables
CandPoint=repmat(Data.xbest,NCandidates,1);
for ii = 1:NCandidates %for each canidate point
    rp=rand(1,Data.dim);
    pertID=zeros(1,Data.dim);
    pertID(rp<P)=1; %find the variables of xbest to perturb
    if sum(pertID)==0 %no variables are to be perturbed -> select one variable at random
        pID=randperm(Data.dim);
        pertID(pID(1))=1;
    end
    for jj=1:Data.dim %for each integer variables
        if  pertID(jj)==1 %check if perturbation wanted
            r=randperm(length(stdev_int));  %select random perturbation range
            p=randn(1);
            sign_p=sign(p);
            if sign_p <0
                perturbation= max(-(Data.xup(jj)-Data.xlow(jj)),...
                    min(-1,-abs(round(stdev_int(r(1))*randn(1)))));
            else
                perturbation= min(Data.xup(jj)-Data.xlow(jj),...
                    max(1,abs(round(stdev_int(r(1))*randn(1)))));
            end
            CandPoint(ii,jj) = max(Data.xlow(jj),...
                min(CandPoint(ii,jj)+perturbation,...
                Data.xup(jj)));
        end
    end
end

%if global sampling desired, add points that are randomly selected fromt the
%whole variable domain
if strcmp(Data.sampling_technique, 'CANDglob')
    %uniformly selected integer points 
    CandPointU= round(repmat(Data.xlow,NCandidates,1) + ...
        rand(NCandidates,Data.dim).*repmat(Data.xup-Data.xlow,NCandidates,1));
    CandPoint=[CandPoint;CandPointU];
end

end% function