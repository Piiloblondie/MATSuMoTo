function CandPoint = Perturbation_SOMI(Data,NCandidates,stdev_int,...
    stdev_cont,P)
%PERTURBATION_SOMI.m generates candidate points for the next sample site by 
%perturbing vaiables of best point found so far. This routine is for
%MIXED INTEGER problems
%Three or four groups are generated: 
%a) only continuous perturbation (small, medium, large)
%b) only integer perturbation (small, medium, large)
%c) continuous and integer perturbation (small, medium, large)
%d) randomly generated points over whole domain (if 'CANDglob' sampling)
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
%Data - structure, containing all problem information
%NCandidates - integer, number of candidate points in each group
%stdev_int - vector, contains perturbation ranges for integers
%stdev_cont - vector, contains perturbation ranges for continuous variables
%P - scalar, perturbation probability
%
%Output:
%CandPoint - matrix with candidates for next function evaluation
%--------------------------------------------------------------------------

%Group1: perturb only continuous variables  
d_cont=length(Data.continuous);   %dimension of continuous variables
%set candidate points same as xbest, and then perturb only continuous
%variables
CandPoint_G1=repmat(Data.xbest,NCandidates,1);
for ii = 1:NCandidates %for each canidate point
    rp=rand(1,d_cont);
    pertID=zeros(1,d_cont);
    pertID(rp<P)=1; %find the variables of xbest to perturb
    if sum(pertID)==0 %no variables are to be perturbed -> select one variable at random
        pID=randperm(d_cont);
        pertID(pID(1))=1;
    end 
    for jj=1:d_cont %for each continuous variables
        if pertID(jj)==1 %check if perturbation wanted
            r=randperm(length(stdev_cont)); %selecet random perturbation range
            CandPoint_G1(ii,Data.continuous(jj)) = max(Data.xlow(Data.continuous(jj)),...
                min(CandPoint_G1(ii,Data.continuous(jj))+stdev_cont(r(1))*randn(1),...
                Data.xup(Data.continuous(jj))));
        end
    end
end


%Group2: perturb only integer variables     
d_int=length(Data.integer);   %dimension of integer variables
%set candidate points same as xbest, and then perturb only integer
%variables
CandPoint_G2=repmat(Data.xbest,NCandidates,1);
for ii = 1:NCandidates %for each canidate point
    rp=rand(1,d_int);
    pertID=zeros(1,d_int);
    pertID(rp<P)=1; %find the variables of xbest to perturb
    if sum(pertID)==0 %no variables are to be perturbed -> select one variable at random
        pID=randperm(d_int);
        pertID(pID(1))=1;
    end 
    
    for jj=1:d_int %for each integer variable
        if pertID(jj)==1 %check if perturbation wanted
            r=randperm(length(stdev_int));  %select random perturbation range
            p=randn(1);
            sign_p=sign(p);
            if sign_p <0
                perturbation= min(-1,-abs(round(stdev_int(r(1))*randn(1))));
            else
                perturbation= max(1,abs(round(stdev_int(r(1))*randn(1))));
            end
            CandPoint_G2(ii,Data.integer(jj)) = max(Data.xlow(Data.integer(jj)),...
                min(CandPoint_G2(ii,Data.integer(jj))+perturbation,...
                Data.xup(Data.integer(jj))));
        end
    end
end

%Group 3: perturb integer and continuous variables
%set candidate points same as xbest, and then perturb integer & continuous
%variables
CandPoint_G3=repmat(Data.xbest,NCandidates,1);
for ii = 1:NCandidates %for each canidate point
    rp=rand(1,Data.dim);
    pertID=zeros(1,Data.dim);
    pertID(rp<P)=1; %find the variables of xbest to perturb
    if sum(pertID)==0 %no variables are to be perturbed -> select one variable at random
        pID=randperm(Data.dim);
        pertID(pID(1))=1;
    end 
    for jj=1:Data.dim %for each variable
        if pertID(jj)==1 && ismember(jj,Data.integer) %integer perturbation
            r=randperm(length(stdev_int)); %selecet random perturbation range
            p=randn(1);
            sign_p=sign(p);
            if sign_p <0
                perturbation= min(-1,-abs(round(stdev_int(r(1))*randn(1))));
            else
                perturbation= max(1,abs(round(stdev_int(r(1))*randn(1))));
            end
            CandPoint_G3(ii,jj) = max(Data.xlow(jj),...
                min(CandPoint_G3(ii,jj)+perturbation,Data.xup(jj)));
        elseif pertID(jj)==1 %continuous perturbation
            r=randperm(length(stdev_cont)); %selecet random perturbation range
            CandPoint_G3(ii,jj) = max(Data.xlow(jj),...
                min(CandPoint_G3(ii,jj)+stdev_cont(r(1))*randn(1),...
                Data.xup(jj)));
        end
    end
end

if strcmp(Data.sampling_technique,'CANDglob')
    %Group4: uniformly sampled points over variable domain
    CandPoint_G4=repmat(Data.xlow,NCandidates,1) + rand(NCandidates,length(Data.xlow)).*...
        repmat(Data.xup-Data.xlow,NCandidates,1);
    %round integer variables
    CandPoint_G4(:,Data.integer)=round(CandPoint_G4(:,Data.integer));
    %put all candidate points together
    CandPoint=[CandPoint_G1;CandPoint_G2;CandPoint_G3; CandPoint_G4];
    %delete unnecessary data
    clear CandPoint_G1 CandPoint_G2 CandPoint_G3 CandPoint_G4 
else
    %put all candidate points together
    CandPoint=[CandPoint_G1;CandPoint_G2;CandPoint_G3];
    %delete unnecessary data
    clear CandPoint_G1 CandPoint_G2 CandPoint_G3
end    

end %function