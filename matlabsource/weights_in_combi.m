function w=weights_in_combi(Intervals,newmodel,Data)
%WEIGHTS_IN_COMBI determines the weights associated with each surrogate
%model in the mixture
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
%Intervals - matrix, contains mixture models with plausibilities, believes,
%and pignistic probabilities
%newmodel - vector, contains indices of the models in the mixture
%Data - structure, contains all optimization problem information
%
%output:
%w - vector with weights of each model contributing to the mixture 
%
%--------------------------------------------------------------------------

I = NaN*ones(length(newmodel),Data.numberOfModels);
for ii=1:length(newmodel)
    I(ii,newmodel(ii))=newmodel(ii);
end

%determine weights of each model in the combined model based on
%pignistic probabilities
w=zeros(length(newmodel),1);
for ii = 1:length(newmodel)
    for jj = 3:size(Intervals,1)
        if  (sum(isnan(Intervals(jj,1:Data.numberOfModels)))==sum(isnan(I(ii,:)))) &&...
                (all(~isnan(Intervals(jj,1:Data.numberOfModels))-~isnan(I(ii,:))==0)) 
            w(ii,1)=Intervals(jj,end);
        end
    end
end
end %function