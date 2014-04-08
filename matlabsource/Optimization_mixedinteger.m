function Solution= Optimization_mixedinteger(Data)
%OPTIMIZATION_MIXEDINTEGER.m is the optimization phase for mixed-integer 
%unconstrained problems: 
%iteratively select new sample site(s), do the expensive
%function evaluation, and update the chosen surrogate model; iterate until
%stopping criterion met
%The algorithm may restart from scratch if it determined being trapped in a
%local optimum
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
%Data - structure with all problem information 
%
%Output:
%Solution - structure containing all information about the solution
%--------------------------------------------------------------------------

%initialize structure that contains all problem-relevant data
Solution.S=[]; %sample points
Solution.Y=[]; %function values corresponding to sample points
Solution.fbest=[]; %best function value found in each trial
Solution.xbest=[]; %best point found in each trial
Solution.fevaltime=[]; %time needed for function evaluation
Solution.trial=0; %number of trials (within maxeval allowed evaluations)
Solution.EvalsEachTrial = []; %number of function evaluations in each trial
UserGivenNumberNewSamples=Data.NumberNewSamples; %the number of desired function evaluations in each iteration

while size(Solution.S,1) < Data.maxeval %do until the budget of allowed function evaluations is exhausted
    %generate initial experimental design
    rankS=0;
    while rankS < Data.dim+1 
        startPoint_cur=Data.starting_point; %user given starting points to add to the design
        Data.S=StartingDesign(Data);  %create initial experimental design
        if isfield (Data,'starting_point') %user gives one or more starting points to add to the initital design
            %check if any given starting point is already constined in initial design
            ii = 1;
            while ii <= size(startPoint_cur,1)
                if any(sum((repmat(startPoint_cur(ii,:),size(Data.S,1),1)-Data.S).^2,2)==0)
                    startPoint_cur(ii,:)=[];
                else
                    ii=ii+1;
                end
            end
            %add remaining starting points to initial experimental design
            Data.S=[startPoint_cur; Data.S];
        end
        %round integer variables
        Data.S(:,Data.integer)=round(Data.S(:,Data.integer));
        rankS = rank([ones(size(Data.S,1)), Data.S]); %compute matrix rank
    end

    %do expensive function evaluations at points in initial experimental design
    if ~Data.parallel_eval
%%---------sequential evaluation-------------------------------------------
        Data.Y=zeros(size(Data.S,1),1); %initialize vector with function values
        for ii = 1:size(Data.S,1)
            fevalst=tic; %record function evaluation time 
            Data.Y(ii,1)=feval(Data.objfunction,Data.S(ii,:)); %expensive evaluation
            Data.fevaltime(ii,1)=toc(fevalst);
        end 
%%---------end of sequential evaluation------------------------------------
    else
%%-------parallel version of doing function evalutaions--------------------
        fvals=zeros(size(Data.S,1),1); %initialize vector for function values
        object=Data.objfunction; %objective function handle
        Samples=Data.S; %points to be evaluated
        Tfeval=zeros(size(Data.S,1),1); %initialize vector with evaluation times
        parfor ii = 1:size(Data.S,1)
            fevalst=tic; %record function evaluation time 
            fvals(ii,1)=feval(object,Samples(ii,:)); %expensive evaluation
            Tfeval(ii,1)=toc(fevalst);
        end  
        Data.Y=fvals;
        Data.fevaltime=Tfeval;
%%--------end of parallel version------------------------------------------
    end

    %determine perturbation probability for each variable when generating
    %candidate points
    if Data.dim>5 %if problem dimension is greater than 5, perturb each variable with probability 5/dimension
        P=max(0.1,5/Data.dim);
    else %if problem dimension is <= 5, perturb all variables with probability 1
        P=1;
    end

    [Data.fbest,idx]=min(Data.Y);  %best function value so far
    Data.xbest=Data.S(idx,:); %best point found so far
    Data.Ymed=Data.Y; %vector where large function values are replaced by median
    MedY=median(Data.Y); %median of function values 
    Data.Ymed(Data.Y>MedY)=MedY; %replace large function values by median, used in calculation of surrogate model parameters

    %compute parameters of (mixture) surrogate models
    [lambda,gamma,mmodel,beta, w_m] = FitSurrogateModel(Data);
   
    %initialize algorithm parameters
    xrange = Data.xup-Data.xlow; %range between lower and upper variable bounds
    minxrange = min(xrange); %smallest variable range
    Data.tolerance = 1e-3; %distance tolerance for 2 points to be considered the same
    stdev_cont(1) = 0.2*minxrange;   %rates for perturbing continuous variables
    stdev_cont(2) = 0.1*minxrange;  %rates for perturbing continuous variables          
    stdev_cont(3) = 0.05*minxrange; %rates for perturbing continuous variables
    stdev_int(1)=max(1, round(stdev_cont(1))); %rates for perturbing integer variables
    stdev_int(2)=max(1, round(stdev_cont(2))); %rates for perturbing integer variables    
    stdev_int(3)=max(1, round(stdev_cont(3))); %rates for perturbing integer variables    
    NCand=125*Data.dim;       %number of candidate points in each group    
    maxshrinkparam = 5;% maximal number of reduction of standard deviation for normal distribution when generating candidate points
    failtolerance = max(5,Data.dim); %max. number of allowed consecutive fails before reducing perturbation range
    succtolerance =3; %max. number of allowed consecutive successes before increasing perturbation range
    shrinkctr = 0;  % number of times perturbation range was decreased
    failctr = 0; % number of consecutive unsuccessful iterations
    succctr=0; % number of consecutive successful iterations
    Data.localminflag=0; %initialize indicator for local minimum
    %iteration loop
    while size(Data.S,1) < (Data.maxeval-size(Solution.S,1)) && ~Data.localminflag
        if strcmp(Data.sampling_technique(1:4),'CAND') %candidate point sampling
            %initialize the number of candidate point groups
            if strcmp(Data.sampling_technique,'CANDglob') % with randomly selected points
                MaxGroupNr=4; %4 groups of candidate points at most
                GroupNr=5;
                RScrit=1.25*ones(1,4);
            elseif strcmp(Data.sampling_technique,'CANDloc') %only local perturbations
                MaxGroupNr=3; %3 groups of candidate points at most
                GroupNr=4;
                RScrit=1.25*ones(1,3);
            end
            Data.NumberNewSamples=UserGivenNumberNewSamples; %reset the number of points to select in each iteration
            Data.NumberNewSamples=min(Data.NumberNewSamples,Data.maxeval-(size(Solution.S,1)+ size(Data.S,1)));
            %determine how many points from each group should be selected
            SampleGroup=zeros(1,Data.NumberNewSamples);          
            for ii = 1:Data.NumberNewSamples
                GroupNr=GroupNr-1;            
                if GroupNr <= 0 && strcmp(Data.sampling_technique,'CANDglob')
                    GroupNr = 4;
                elseif GroupNr <= 0 
                    GroupNr = 3;            
                end
                SampleGroup(ii)=GroupNr;
            end

            %create candidate points by perturbing discrete, continuous, discrete
            %and continuous variables (and by randomly sampling points
            %over domain)
            CandPoint = Perturbation_SOMI(Data,NCand,stdev_int,stdev_cont,P);
            %cycling weights for each group - determine response surface
            %weight for each sample point
            xselected=[];
            weights={};
            %determine how many points to be selected from each group
            %and what the weights for the criteria should be
            for ii = 1:MaxGroupNr
                Data.NumberNewSamples=length(find(SampleGroup==ii)); %number of points to be selected from group ii
                for jj = 1:Data.NumberNewSamples
                    RScrit(ii)=RScrit(ii)-.25;
                    if RScrit(ii) < 0
                        RScrit(ii) = 1;
                    end
                    weights{ii}(jj) =RScrit(ii);
                end
                if Data.NumberNewSamples > 0
                    %select new points to evaluate from each group
                    xselected_add=SamplePointSelection(Data,CandPoint((ii-1)*NCand+1:ii*NCand,:), weights{ii},...
                        lambda,gamma,mmodel,beta, w_m);
                else
                    xselected_add=[];
                end
                xselected=[xselected;xselected_add];
            end
            clear CandPoint
        else %minimum of response surface sampling
            xselected = SurfMin(Data,lambda,gamma, mmodel,beta, w_m);
        end
    
        no_xnew=size(xselected,1); %number of new sample sites   
        % perform expensive function evaluation at the selected point
        Fselected=zeros(no_xnew,1);
        if no_xnew >1 && Data.parallel_eval
%%-------parallel version of doing function evalutaions--------------------
              object=Data.objfunction; %objective function handle
              Samples=xselected; %points to be evaluated
              Tfeval=zeros(no_xnew,1); %initialize vector with function evaluation times
              parfor ii = 1:no_xnew %parallel for
                  fevalst=tic; %record function evaluation time 
                  Fselected(ii,1)=feval(object,Samples(ii,:)); %expensive evaluation
                  Tfeval(ii,1)=toc(fevalst);
              end  
              Data.fevaltime=[Data.fevaltime;Tfeval];
%%--------end of parallel version------------------------------------------
        else
%%---------------sequentially evaluating new points------------------------
            for ii=1:no_xnew
                fevalst=tic; %record evaluation time
                Fselected(ii,1) = feval(Data.objfunction,xselected(ii,:)); %expensive evaluation
                Data.fevaltime=[Data.fevaltime;toc(fevalst)];
            end
%%-----------end of sequential evaluation----------------------------------
        end

        %check if new points are improvements
        Success=0;
        for ii = 1:no_xnew
            if Fselected(ii) < Data.fbest %point has better function value than best point found so far
                if (Data.fbest - Fselected(ii)) > (1e-3)*abs(Data.fbest) %improvement must be more than 0.1% 
                    Success=true; %better point found
                end
                Data.xbest = xselected(ii,:); %update best point found so far
                Data.fbest = Fselected(ii); %update best function value found so far
            end
        end

        if Success %at least one new point is better than best point found so far
            failctr = 0;  %re-initialize fail-counter
            succctr=succctr+1; %update success counter
            shrinkctr = 0; %re-initialize shrink-counter
        else 
            failctr = failctr + 1;  %update fail-counter
            succctr=0; %re-initialize success counter
        end

        % check if algorithm is in a local minimum
        shrinkflag = 1;      
        if failctr >= failtolerance %check how many consecutive failed improvement trials
            if shrinkctr >= maxshrinkparam
                shrinkflag = 0;  %re-initialize shrink indicator
                disp('Stopped reducing perturbation radius because the maximum number of reductions has been reached.');
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
                stdev_cont_new = max(stdev_cont/2,Data.tolerance);
                if all(stdev_cont_new ==stdev_cont)
                    disp('Cannot reduce continuous perturbation any further!');
                else
                    stdev_cont=stdev_cont_new;
                    disp('Reducing continuous perturbation by half!');
                end
            else  %local min detected
                Data.localminflag = 1;
                disp('Algorithm is probably in a local minimum! Restarting the algorithm from scratch.');
            end
        end

        if succctr>=succtolerance %check if number of consecutive improvements is large enough        
            stdev_int=round(min(stdev_int*2, max(Data.xup(Data.integer)-Data.xlow(Data.integer))));%increase search radius integer
            stdev_cont=min(stdev_cont*2, max(Data.xup(Data.continuous)-Data.xlow(Data.continuous)));%increase search radius continuous
            succctr=0;  %re-initialize success counter
            shrinkctr = shrinkctr - 1; %update shrink counter
        end
    
        Data.S = [Data.S; xselected]; %update matrix of sample sites
        Data.Y = [Data.Y; Fselected]; %update vector of objective function values
        Data.Ymed=Data.Y; %vector where large function values are replaced by median
        MedY=median(Data.Y); %median of function values 
        Data.Ymed(Data.Y>MedY)=MedY; %replace large function values by median, used in calculation of surrogate model parameters     
        
        if size(Solution.S,1)+ size(Data.S,1)< Data.maxeval && Data.localminflag == 0
            %update model parameters
            [lambda,gamma,mmodel,beta, w_m] = FitSurrogateModel(Data);
            fprintf('Total number of function evaluations: %4.0f;\n',size(Solution.S,1) +  size(Data.S,1))
            fprintf('Best function value in this trial: %f\n',Data.fbest)
            if ~isempty(Solution.fbest) 
                fprintf('Best function value over all previous trials: %f\n',min(Solution.fbest))
            else
                fprintf('Best function value over all previous trials: %f\n',nan)
            end
            save('Results','Data','Solution'); %save intermediate results
        elseif Data.localminflag == 1 && size(Solution.S,1) + size(Data.S,1) +...
                Data.number_startpoints +size(Data.starting_point,1)>= Data.maxeval 
            %local min encountered and not sufficiently many points in
            %budget to restart from scratch; do remaining steps as if no
            %local min detected
            if size(Solution.S,1) + size(Data.S,1)  >= Data.maxeval
                % budget of function evaluations exhausted
                Solution.S=[Solution.S;Data.S]; %update sample site matrix
                Solution.Y=[Solution.Y; Data.Y]; %update function value vector
                Solution.fbest=[Solution.fbest; Data.fbest];  %update best function value found
                Solution.xbest=[Solution.xbest; Data.xbest];  %update best point found
                Solution.fevaltime=[Solution.fevaltime; Data.fevaltime]; %update function evaluation time vector
                Solution.trial = Solution.trial + 1; %update number of restarts from scratch
                Solution.EvalsEachTrial = [Solution.EvalsEachTrial;length(Data.Y)]; %update number of evaluations in each trial
                Data.Ymed=[];  %re-initialize median function value vector
                Data.S=[]; %re-initialize sample site matrix 
                Data.Y=[]; %re-initialize function value vector
                Data.fevaltime=[]; %re-initialize function evaluation time vector
                Data.xbest=[]; %re-initialize best point found so far
                fprintf('Total number of function evaluations: %4.0f;\n',size(Solution.S,1) +  size(Data.S,1))
                fprintf('Best function value in this trial: %f\n',Data.fbest)
                Data.fbest=[];%re-initialize best function value found so far
            else           
                Data.localminflag =0; %re-initialize local minimum indicator (there are not enough evaluations left for a restart)
                [lambda,gamma,mmodel,beta, w_m] = FitSurrogateModel(Data); %update surrogate model parameters
                fprintf('Total number of function evaluations: %4.0f;\n',size(Solution.S,1) +  size(Data.S,1))
                fprintf('Best function value in this trial: %f\n',Data.fbest)
            end
            if ~isempty(Solution.fbest) 
                fprintf('Best function value over all previous trials: %f\n',min(Solution.fbest))
            else
                fprintf('Best function value over all previous trials: %f\n',nan)
            end
            save('Results','Data','Solution');%save intermediate results            
        elseif Data.localminflag == 1 || size(Solution.S,1) + size(Data.S,1) >= Data.maxeval
            %local min detected or budget of function evaluations exhausted
            Solution.S=[Solution.S;Data.S]; %update sample site matrix
            Solution.Y=[Solution.Y; Data.Y]; %update function value vector
            Solution.fbest=[Solution.fbest; Data.fbest]; %update best function value found
            Solution.xbest=[Solution.xbest; Data.xbest]; %update best point found
            Solution.fevaltime=[Solution.fevaltime; Data.fevaltime]; %update function evaluation time vector
            Solution.trial = Solution.trial + 1; %update number of restarts from scratch
            Solution.EvalsEachTrial = [Solution.EvalsEachTrial;length(Data.Y)]; %update number of evaluations in each trial
            Data.Ymed=[];  %re-initialize median function value vector
            Data.S=[]; %re-initialize sample site matrix 
            Data.Y=[]; %re-initialize function value vector
            Data.fevaltime=[]; %re-initialize function evaluation time vector            
            Data.xbest=[]; %re-initialize best point found so far
            fprintf('Total number of function evaluations: %4.0f;\n',size(Solution.S,1) +  size(Data.S,1))
            fprintf('Best function value in this trial: %f\n',Data.fbest)
            Data.fbest=[]; %re-initialize best function value found so far
            if ~isempty(Solution.fbest) 
                fprintf('Best function value over all previous trials: %f\n',min(Solution.fbest))
            else
                fprintf('Best function value over all previous trials: %f\n',nan)
            end
            save('Results','Data','Solution'); %save intermediate results          
        end
    end
end
end %function 