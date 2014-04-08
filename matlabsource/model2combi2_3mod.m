function Ymixtest=model2combi2_3mod(Intervals,yPred_all)
%MODEL2COMBI2_3MOD.m computes weights for models in 3-model mixtures and 
%the predictions of the mixture model
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
%Intervals - matrix with believes, plausibilities, pignistic probabilities
%            for each model in the mixture
%yPred_all - matrix with predictions of single surrogate models
%
%output:
%Ymixtest - mixture model prediction
%--------------------------------------------------------------------------

a1=find(Intervals(:,1)==1);
a2=find(Intervals(:,2)==2);
a3=find(Intervals(:,3)==3);

Inew=Intervals([a1,a2,a3],[1:3,end]);

Ymixtest=zeros(size(yPred_all,1),4); %(m x 11) matrix combinations in colns: 1/2 1/3 1/4 2/3 2/4 3/4 1/2/3 1/2/4 1/3/4 2/3/4 1/2/3/4
% contains predictions of mixed models at sampled sites
Mweight2=[];
for ii = 1:2
    for jj=ii+1:3
        weights=Inew([ii,jj],end)./sum(Inew([ii,jj],end));
        Mweight2=[Mweight2; ii jj weights(:)'];    %weights for combining 2 models   
    end
end
for ii = 1:size(Mweight2,1)
    Ymixtest(:,ii)=Mweight2(ii,3)*yPred_all(:,Mweight2(ii,1))+Mweight2(ii,4)*yPred_all(:,Mweight2(ii,2)); %predictions at sample points in LOOCV
end

%mixture model prediction
Mweight3=Inew(:,end)'./sum(Inew(:,end));
Ymixtest(:,end)=Mweight3(1,1)*yPred_all(:,1)+Mweight3(1,2)*yPred_all(:,2)+Mweight3(1,3)*yPred_all(:,3);

end %function