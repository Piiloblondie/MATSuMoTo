function xselected=SamplePointSelection(Data,CandPoint,RScrit,lambda,...
    gamma,mmodel,beta, w_m)
%SAMPLEPOINTSELECTION.m is used to select new evaluation points from a set
%of candidate points. This is for candidate point sampling only. A scoring 
%of candidate points is used: Compute the distance of all candidates to 
%already evaluated points and predict the objective function values for
%each candidate point with the surrogate model. Both criteria (distance and
%predicted objective function value are scaled to [0,1], and the candidate 
%with the best weighted score of both criteria becomes the new sample 
%point. If there are more than one new sample point to be selected, the 
%distances of the candidate points to the already selected candidate points
%have to be taken into account.
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
%Data: structure with all problem information
%CandPoint: matrix with candidate points
%RScrit: weight for response surface criterion
%lambda,gamma: parameters of RBF model ([] if RBF not used)
%mmodel: MARS model parameters (from ARESlab toolbox, [] if not used)
%beta: parameters of polynomial regression model ([] if not used)
%w_m - weights for models in mixtures 
%
% Output:
% xselected: matrix with all seleced points for next evaluation
%--------------------------------------------------------------------------  

% predict objective function value at the candidate points
CandValue=PredictFunctionValues(Data,CandPoint,lambda,gamma,beta,mmodel,w_m);
%compute the distance between candidate points and all already sampled points
NormValue=pdist2(CandPoint,Data.S)';  %pdist2 is MATLAB built-in function

%scale predicted objective function values to [0,1]
MinCandValue = min(CandValue); %find min of predIcted objective function value
MaxCandValue = max(CandValue); %find maximum of predicted objective function value
if MinCandValue == MaxCandValue  %compute scaled objective function value scores
    ScaledCandValue=ones(length(CandValue),1);
else
    ScaledCandValue = (CandValue-MinCandValue)/(MaxCandValue-MinCandValue);
end

if Data.NumberNewSamples == 1 %only one point is selected
    valueweight=RScrit(1); %weight for response surface criterion
    %scale distances to already evaluated points to [0,1]
    CandMinDist = (min(NormValue,[],1))'; %minimum distance of every candidate point to already sampled points  
    MaxCandMinDist = max(CandMinDist); %maximum of distances
    MinCandMinDist = min(CandMinDist); %minimum of distances
    if  MaxCandMinDist ==MinCandMinDist  %compute distance criterion scores
        ScaledCandMinDist =ones(length(CandMinDist),1);
    else
        ScaledCandMinDist = (MaxCandMinDist-CandMinDist)/(MaxCandMinDist-MinCandMinDist);
    end

    %compute weighted score for all candidates
    CandTotalValue = valueweight*ScaledCandValue + (1 - valueweight)*ScaledCandMinDist;

    %assign inf scores to candidate points that are too close to already sampled
    %points
    CandTotalValue(CandMinDist < Data.tolerance) = Inf; 

    %find candidate with best score -> becomes new sample point
    [MinCandTotalValue,selindex] = min(CandTotalValue);
    if MinCandTotalValue == inf %all candidates are too close to already evaluated points
        xselected=[];
    else
        xselected(1,:) = CandPoint(selindex,:);  %new sample point
    end

else  %more than one new sample point wanted
    xselected=[]; %initialize matrix for new evaluation points
    NrSelected=1;
    for ii =1:Data.NumberNewSamples 
        valueweight=RScrit(ii); %response surface criterion weight        
        if ii == 1 %select first candidate point
            %scale distance values to [0,1]
            CandMinDist = (min(NormValue,[],1))'; %minimum distance of every candidate point to already sampled points  
            MaxCandMinDist = max(CandMinDist); %maximum of distances
            MinCandMinDist = min(CandMinDist); %minimum of distances
            if  MaxCandMinDist ==MinCandMinDist  %compute distance criterion scores
                ScaledCandMinDist =ones(length(CandMinDist),1);
            else
                ScaledCandMinDist = (MaxCandMinDist-CandMinDist)/(MaxCandMinDist-MinCandMinDist);
            end

            %compute weighted score for all candidates
            CandTotalValue = valueweight*ScaledCandValue + (1 - valueweight)*ScaledCandMinDist;

            %assign inf scores to candidate points that are too close to already sampled
            %points
            CandTotalValue(CandMinDist < Data.tolerance) = Inf; 
            if all(CandTotalValue==inf)
                xselected=[];
                return %all candidate points have already been sampled. 
            end
            %find candidate with best score -> becomes new sample point
            [MinCandTotalValue,selindex] = min(CandTotalValue);

            if MinCandTotalValue == inf
                xselected=[xselected;[]];
            else
                xselected = [xselected;CandPoint(selindex,:)];  %new sample point
                NrSelected=NrSelected+1;
            end
        else    
            %compute distance of all candidate points to the previously selected
            %candidate point
            if NrSelected > 1
                NormValueP=sqrt(sum((repmat(xselected(NrSelected-1,:),size(CandPoint,1),1)-CandPoint).^2,2));
                NormValue=[NormValue;NormValueP']; %augment distance matrix
            end

            %re-scale distance values to [0,1]
            CandMinDist = (min(NormValue,[],1))'; %minimum distance of every candidate point to already sampled points  
            MaxCandMinDist = max(CandMinDist); %maximum of distances
            MinCandMinDist = min(CandMinDist); %minimum of distances
            if  MaxCandMinDist ==MinCandMinDist  %compute distance criterion scores
                ScaledCandMinDist =ones(length(CandMinDist),1);
            else
                ScaledCandMinDist = (MaxCandMinDist-CandMinDist)/(MaxCandMinDist-MinCandMinDist);
            end

            %compute weighted score for all candidates
            CandTotalValue = valueweight*ScaledCandValue + (1 - valueweight)*ScaledCandMinDist;    
            %assign inf values to points that are too close to already
            %evaluated/chosen points
            CandTotalValue(CandMinDist < Data.tolerance)=Inf;

            [MinCandTotalValue,selindex] = min(CandTotalValue); %find best candidate
            if MinCandTotalValue == inf
                xselected=[xselected;[]];
            else
                xselected = [xselected;CandPoint(selindex,:)];  %new sample point
                NrSelected=NrSelected+1;
            end
        end
    end
end
ii=1; %delete identical selected sample points
while ii <= size(xselected,1)-1 
    jj=ii+1;
    while jj <= size(xselected,1)
        if sqrt(sum((xselected(ii,:)-xselected(jj,:)).^2)) < Data.tolerance
            xselected(jj,:)=[];
        else
            jj=jj+1;
        end
    end
    ii=ii+1;
end

if isempty(xselected) %regenerate random sample points if no new point selected
    while 1
        xselected = Data.xlow+(Data.xup-Data.xlow).*rand(1,Data.dim);
        xselected(Data.integer)=round(xselected(Data.integer));
        dist_new=min(pdist2(xselected,Data.S));
        if dist_new>Data.tolerance
            break
        end
    end
end    
end%function