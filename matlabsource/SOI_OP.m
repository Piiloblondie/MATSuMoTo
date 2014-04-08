function Data = SOI_OP(Data,P,stdev_int,Solution)
%SOI_OP.m is the optimization phase for purely integer problems. We
%iteratively select new sample sites for doing expensive evaluations and
%update the response surface until the budget of allowable evaluations is
%used up or we are in a local minimum in which case we start from scratch.
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
% Data - structure, contains all problem information
% P - perturbation probability for creating canididate points
% stdev_int - perturbation range for integer variables
% Solution - structure with data from previous trials
%
%Output:
%Data - updated structure with problem information
%--------------------------------------------------------------------------
[Data.fbest,idx]=min(Data.Y); %find best function value
Data.xbest=Data.S(idx,:);     %find best point

Data.Ymed=Data.Y; %vector where large function values are replaced by median
MedY=median(Data.Y); %median of function values 
Data.Ymed(Data.Y>MedY)=MedY; %replace large function values by median, used in calculation of surrogate model parameters

%algorithm parameters
NCandidates=500*Data.dim;   %number of candidate points in each group 
RScrit=1.25; %initialize weight for response surface criterion
%compute parameters of response surface
[lambda,gamma,mmodel,beta, w_m] = FitSurrogateModel(Data);
   
%initialize parameters for determining whether or not local minimum has
%been found
maxshrinkparam = 5;% maximal number of reductions of standard deviation for normal distribution when generating candidate points
failtolerance = max(5,Data.dim); %max. number of allowed consecutive fails before reducing perturbation range
succtolerance =3; %max. number of allowed consecutive successes before increasing perturbation range
shrinkctr = 0; % number of times perturbation range was decreased
failctr = 0; % number of consecutive unsuccessful iterations
succctr=0; % number of consecutive successful iterations
Data.localminflag=0; %initialize indicator for local minimum

% iterate until max number of function evaluations reached or local min
% found
while size(Solution.S,1)+ size(Data.S,1) < Data.maxeval && ~Data.localminflag 
    %set the weights for the response surface criterion for each sample
    %point
    if strcmp(Data.sampling_technique(1:4), 'CAND') %randomized smapling startegy
        weights=zeros(1,Data.NumberNewSamples); %initialize weight vector for criteria
        Data.NumberNewSamples=min(Data.NumberNewSamples,Data.maxeval-(size(Solution.S,1)+ size(Data.S,1)));
        % determine which weights to use    
        for ii = 1:Data.NumberNewSamples
            RScrit=RScrit-.25;
            if RScrit < 0
                RScrit = 1;
            end
            weights(ii) =RScrit;
        end
        %generate candidate points          
        CandPoint=Perturbation_SOI(Data,NCandidates,stdev_int,P);
        %select new sample points
        xselected=SamplePointSelection(Data,CandPoint, weights,...
            lambda,gamma,mmodel,beta, w_m);
        clear CandPoint 
    else %sample at minimum point of response surface
        xselected = SurfMin(Data,lambda,gamma, mmodel,beta, w_m);
    end
        
    no_xnew=size(xselected,1); %number of new sample sites   
    %  perform expensive function evaluation at the selected points
    Fselected=zeros(no_xnew,1); %initialize vector with new function values
    if no_xnew >1 && Data.parallel_eval
%%-------parallel version of doing function evalutaions--------------------
        object=Data.objfunction; %objective function handle
        Samples=xselected; %points to be evaluated
        Tfeval=zeros(no_xnew,1); %initialize vector with function evaluation times
        parfor ii = 1:no_xnew
            fevalst=tic; %record function evaluation time 
            Fselected(ii,1)=feval(object,Samples(ii,:)); %expensive evaluation
            Tfeval(ii,1)=toc(fevalst);
        end  
        Data.fevaltime=[Data.fevaltime;Tfeval];
%%--------end of parallel version------------------------------------------
    else          
%%-------sequentially evaluating new points--------------------------------
        for ii=1:no_xnew
            fevalst=tic; %record evaluation time
            Fselected(ii,1) = feval(Data.objfunction,xselected(ii,:)); %expensive evaluation
            Data.fevaltime=[Data.fevaltime;toc(fevalst)];
        end
%%--------end of sequential evaluation-------------------------------------
    end

    %check if new points are feasible and improvements
    Success=false;
    for ii = 1:length(Fselected) 
        if Fselected(ii) < Data.fbest %point has better function value then best point found so far
            if (Data.fbest - Fselected(ii)) > (1e-3)*abs(Data.fbest) %improvement must be more than 0.1%
                Success=true;
            end
            Data.xbest = xselected(ii,:); %update best point found so far
            Data.fbest =  Fselected(ii); %update best function value found so far
        end
    end
    
    if Success %better point found in this iteration
        failctr = 0; %re-initialize fail-counter
        succctr=succctr+1; %update success counter
        shrinkctr = 0; %re-initialize shrink-counter
    else
        failctr = failctr+1; %update fail-counter
        succctr=0; %re-initialize success counter
    end

    % check if algorithm is in a local minimum
    shrinkflag = 1;      
    if failctr >= failtolerance %check how many consecutive failed improvement trials
        if shrinkctr >= maxshrinkparam %reduced perturbaltion radius as often as possible
            shrinkflag = 0; %re-initialize shrink indicator
            disp('Stopped reducing perturbation radius because the maximum number of reduction has been reached.');
        end
        failctr = 0; %re-initialize fail-counter
        if shrinkflag == 1 %reduce perturbation range
            shrinkctr = shrinkctr + 1; %update shrink-counter
            stdev_int_new = max(round(stdev_int/2),1);
            if all(stdev_int_new ==stdev_int)
                disp('Cannot reduce integer perturbation any further!');
            else
                stdev_int=stdev_int_new;
                disp('Reducing integer perturbation by half!');
            end
        else  %local minimum
            Data.localminflag = 1;
            disp('Algorithm is probably in a local minimum! Restarting the algorithm from scratch.');
        end
    end
    
    if succctr>=succtolerance %check if number of consecutive improvements is large enough        
        stdev_int=round(min(stdev_int*2, max(Data.xup-Data.xlow)));%increase search radius
        succctr=0; %re-initialize success counter
        shrinkctr = 0; %re-initialize shrink counter
    end
    
    Data.S = [Data.S; xselected]; %update sample site matrix
    Data.Y = [Data.Y; Fselected]; %update function value vector
    Data.Ymed=Data.Y; %vector where large function values are replaced by median
    MedY=median(Data.Y); %median of function values 
    Data.Ymed(Data.Y>MedY)=MedY; %replace large function values by median, used in calculation of surrogate model parameters

    %update model parameters
    if size(Solution.S,1)+size(Data.S,1)< Data.maxeval && Data.localminflag == 0
        [lambda,gamma,mmodel,beta, w_m] = FitSurrogateModel(Data);
    end
    if Data.localminflag == 1 && size(Solution.S,1) + size(Data.S,1) +...
        Data.number_startpoints+size(Data.starting_point,1) >= Data.maxeval
        %local minimum encountered, but budget of ramaining evaluations is
        %not sufficient to restart
        Data.localminflag=0;
        [lambda,gamma,mmodel,beta, w_m] = FitSurrogateModel(Data); %update model parameters
    end
    
    fprintf('Total number of function evaluations: %4.0f;\n',size(Solution.S,1) +  size(Data.S,1))
    fprintf('Best function value in this trial: %f\n',Data.fbest)
    if ~isempty(Solution.fbest) && ~all(isnan(Solution.fbest))
        fprintf('Best function value over all trials: %f\n',min(Solution.fbest))
    elseif isempty(Solution.fbest) 
        fprintf('Best function value over all previous trials: %f\n',nan)
    end
    save('Results','Data','Solution');%save intermediate results
end

end % function