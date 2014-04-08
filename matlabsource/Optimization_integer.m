function Solution= Optimization_integer(Data)
%OPTIMIZATION_INTEGER.m is the optimization phase for pure integer 
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
%Data - structure, contains all information about the optimization problem
%Output:
%Solution -  structure with all problem-relevant data
%--------------------------------------------------------------------------

%initialize structure that contains all problem-relevant data
Solution.S=[]; %sample points
Solution.Y=[]; %function values corresponding to sample points
Solution.fbest=[]; %best function value found in each trial
Solution.xbest=[]; %best point found in each trial
Solution.fevaltime=[]; %time needed for function evaluation
Solution.trial=0; %number of trials (within maxeval allowed evaluations)
Solution.EvalsEachTrial = []; %number of function evaluations in each trial

% do until the budget of allowed function evaluations is exhausted
while size(Solution.S,1) < Data.maxeval 
    Data.tolerance=1; %toleranc below which points are identical
    %regenerate initial experimental design until sample site matrix has
    %correct rank (rank[1,S] = dimension+1)
    rankS=0;
    while rankS < Data.dim+1
        startPoint_cur=Data.starting_point; %user supplied starting points
        Data.S=StartingDesign(Data); %create initial experimental design
        if  isfield (Data,'starting_point') %user gives one or more starting points to add to the initital design
            %check if any given starting point is already contained in initial design
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
        Data.S=round(Data.S);     %round integer variables
        rankS = rank([ones(size(Data.S,1)), Data.S]); %compute matrix rank
    end

%%-------------sequential evaluation of points in initial design-----------
    if ~Data.parallel_eval
        Data.Y=zeros(size(Data.S,1),1); %initialize vector with function values
        for ii = 1:size(Data.S,1)
            fevalst=tic; %record function evaluation time 
            Data.Y(ii,1)=feval(Data.objfunction,Data.S(ii,:)); %expensive evaluation
            Data.fevaltime(ii,1)=toc(fevalst);
        end 
%%-----------------end of sequential evaluation----------------------------
    else
%%-----parallel version of doing function evalutaions----------------------
        fvals=zeros(size(Data.S,1),1); %initialize vector for function values
        object=Data.objfunction;  %objective function handle
        Samples=Data.S;  %points to be evaluated in parallel
        Tfeval=zeros(size(Data.S,1),1); %initialize vector for evaluation times
        parfor ii = 1:size(Data.S,1) %parallel for
            fevalst=tic; %record function evaluation time 
            fvals(ii,1)=feval(object,Samples(ii,:)); %expensive function evaluations
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
    stdev_int=[1,2,3]; %perturbation ranges 

    [Data.fbest,idx]=min(Data.Y); %best function value in initial experimental design
    Data.xbest=Data.S(idx,:);  %best point in initial experimental design
    
    %optimization phase 
    if size(Solution.S,1) + size(Data.S,1) < Data.maxeval 
        Data = SOI_OP(Data,P,stdev_int,Solution);
    end
    
    % after algorithm stops, save all data in structure Solution. If all
    % function evaluations are used up, stop, otherwise start new
    % optimization trial
    if size(Solution.S,1) + size(Data.S,1) + Data.number_startpoints +size(Data.starting_point,1) >= Data.maxeval
        Data.number_startpoints = Data.maxeval -(size(Solution.S,1) + size(Data.S,1)) - size(Data.starting_point,1);
    end
    Solution.S=[Solution.S;Data.S]; %update sample site matrix
    Solution.Y=[Solution.Y; Data.Y]; %update function value vector
    Solution.fbest=[Solution.fbest; Data.fbest]; %update best function value found
    Solution.xbest=[Solution.xbest; Data.xbest]; %update best point found
    Solution.fevaltime=[Solution.fevaltime; Data.fevaltime]; %update function evaluation time vector
    Solution.trial = Solution.trial + 1; %update number of restarts from scratch
    Solution.EvalsEachTrial = [Solution.EvalsEachTrial;length(Data.Y)]; %update number of evaluations in each trial
    Data.S=[]; %re-initialize sample site matrix 
    Data.Y=[]; %re-initialize function value vector
    Data.fevaltime=[]; %re-initialize function evaluation time vector
    Data.fbest=[]; %re-initialize best function value found so far
    Data.xbest=[]; %re-initialize best point found so far
end
end %function