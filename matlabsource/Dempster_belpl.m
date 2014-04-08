function IntervalsD=Dempster_belpl(BPA_dempster,Data)
%DEMPSTER_BELPL.m calculates believes, plausibilities, and pignistic 
%probabilities of focal elements based on Dempster's rule of combination.
%See literature on Dempster-Shafer theory for definition of these values.
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
%BPA_dempster - matrix with basic probability assignments
%Data - structure with all model information
%Output:
%IntervalsD - matrix with belive, plausibility and pignistic probability
%values
%--------------------------------------------------------------------------

b=(1:Data.numberOfModels);

%% believes 
belD=zeros(size(BPA_dempster,1),1);
for ii = 1:size(BPA_dempster,1)
    cur1=b(~isnan(BPA_dempster(ii,1:Data.numberOfModels)));
    for jj = 1:size(BPA_dempster,1)
        cur2=b(~isnan(BPA_dempster(jj,1:Data.numberOfModels)));
        if all(ismember(cur2,cur1))
            belD(ii)=belD(ii)+BPA_dempster(jj,end);
        end
    end
end

%% plausibilities
plD=zeros(size(BPA_dempster,1),1);
for ii = 1:size(BPA_dempster,1)
    cur1=b(~isnan(BPA_dempster(ii,1:Data.numberOfModels)));
    for jj = 1:size(BPA_dempster,1)
        cur2=b(~isnan(BPA_dempster(jj,1:Data.numberOfModels)));
        if ~isempty(intersect(cur2,cur1))
            plD(ii)=plD(ii)+BPA_dempster(jj,end);
        end
    end
end
    
%% Pignistic probabilities (see FLOREA JOUSSELME BOSSE GRENIER: robust combination rules for evidence theory)
betP_D=zeros(size(BPA_dempster,1),1);
for ii = 1:size(BPA_dempster,1)
    cur1=b(~isnan(BPA_dempster(ii,1:Data.numberOfModels)));
    for jj = 1:size(BPA_dempster,1)
        cur2=b(~isnan(BPA_dempster(jj,1:Data.numberOfModels)));
        if ~isempty(intersect(cur2,cur1))
            betP_D(ii)=betP_D(ii)+BPA_dempster(jj,end)*length(intersect(cur2,cur1))/length(cur2);
        end
    end
end

%% believe-plausibility Intervals
IntervalsD=[BPA_dempster(:,1:Data.numberOfModels) belD plD betP_D];
IntervalsD=sortrows(IntervalsD,-(Data.numberOfModels+2));

end %function