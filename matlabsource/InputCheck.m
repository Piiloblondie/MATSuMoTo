function Data = InputCheck(data_file,varargin)
%INPUTCHECK.m goes through all provided user input and checks if it is
%correct and/or exists and assigns default values to variables that have
%not been set by the user. A pool of MATLAB workers is opened if the
%Parallel Computing Toolbox is installed.
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
%data_file - string with the name of the file that contains all 
%optimization probleminformation
%varargin - cell, containing all other possibly provided input data from
%the user
%
%Output:
%Data - structure with all information about the problem:
% maxeval - maximum number of allowed function evaluations
% surrogate_model - string with name of surrogate model to be used
% sampling_technique -  string with sampling technique to be used
% initial_design - string with initial experimental design to be used
% number_startpoints - number of points for initial experimental design
% starting_point - matrix with starting points to be added to initial
%                 experimental design
% NumberNewSamples - number of points to be selected in every iteration for
%                   doing expensive function evaluation
%--------------------------------------------------------------------------


%% ---------- START INPUT CHECK-------------------
if nargin < 1 ||  ~ischar(data_file) %check if input data file present
    error(['You have to supply a file name with your data. \n' ...
    'See example files for information how to define problems.'])
else
    Data=feval(data_file); % load problem data
end

%check if optimization problem dimension is defined
if ~isfield(Data,'dim')
    error('You must provide the problem dimension.')
elseif abs(Data.dim-round(Data.dim))>0 || Data.dim <= 0 
    error('Dimension must be positive integer.')
end

%check if all upper and lower bounds are given
if length(Data.xlow) ~= Data.dim || length(Data.xup) ~= Data.dim
    error('Vector length of lower and upper bounds must equal problem dimension')
end

%check if lower bounds < upper bounds
if any(Data.xlow >= Data.xup)
    error('Lower bounds have to be lower than upper bounds.')
end

%check if continuous and integer variables are defined
if ~isfield(Data,'integer')
    error(['You must provide the indices of integer variables. ',...
        'Use Data.integer = [] if the problem has continuous variables only.'])
end
if ~isfield(Data,'continuous')
    error(['You must provide the indices of continuous variables. '...
        'Use Data.continuous = [] if the problem has integer variables only.'])
end

%check how many input arguments are given and use default values  for 
%arguments not provided
if length(varargin)<1 %default settings for everything
    fprintf(['No maximal number of allowed function evaluations defined. \n ',...
             'No surrogate model defined. \n '...
             'No sampling strategy defined. \n '...
             'No initial design defined. \n '...
             'No number of starting points defined. \n '...
             'No starting points defined. \n '...
             'No number of new samples defined. \n '...
             'Using default values.\n']);
    Data.maxeval= 20*Data.dim;  %if maximum number of function evaluations not defined, set to 400
    Data.surrogate_model='RBFcub';   %if no surrogate model defined, use cubic RBF
    Data.sampling_technique='CANDglob';  %if no sampling method definded, use candidate points
    Data.initial_design='LHS';       %if no initial experimental design strategy defined, use Latin hypercube
    Data.number_startpoints=2*(Data.dim+1); %default number of starting points 
    Data.starting_point=[]; %by default no starting point defined
    Data.NumberNewSamples=1;  %by default only one point selected in every iteration
elseif length(varargin) < 2
    fprintf(['No surrogate model defined. \n '...
            'No sampling strategy defined. \n '...
            'No initial design defined. \n '...
            'No number of starting points defined.\n '...
            'No starting points defined. \n '...
            'No number of new samples defined. \n '...
            'Using default values.\n']);
    Data.maxeval=varargin{1};    %user defined maximum number of allowed function evaluations
    Data.surrogate_model='RBFcub';   %if no surrogate model defined, use cubic RBF
    Data.sampling_technique='CANDglob';  %if no sampling method definded, use candidate points
    Data.initial_design='LHS';      %if no initial experimental design strategy defined, use Latin hypercube
    Data.number_startpoints=2*(Data.dim+1); %default number of starting points
    Data.starting_point=[]; %by default no starting point defined
    Data.NumberNewSamples=1;  %by default only one point selected in every iteration
elseif length(varargin) < 3
    fprintf(['No initial design strategy defined. \n '...
            'No sampling strategy defined. \n '...
            'No number of starting points defined. \n '...
            'No starting points defined. \n '...
            'No number of new samples defined. \n '...
            'Using default values.\n']);
    Data.maxeval=varargin{1};    %user defined maximum number of allowed function evaluations
    Data.surrogate_model=varargin{2}; %user defined surrogate model type
    Data.sampling_technique='CANDglob';  %if no sampling method definded, use candidate points
    Data.initial_design='LHS';      %if no initial experimental design strategy defined, use Latin hypercube
    Data.number_startpoints=2*(Data.dim+1); %default number of starting points
    Data.starting_point=[]; %by default no starting point defined
    Data.NumberNewSamples=1;  %by default only one point selected in every iteration
elseif length(varargin) < 4
    fprintf(['No initial design strategy defined. \n ',...
            'No number of starting points defined. \n '...
            'No starting points defined. \n '...
            'No number of new samples defined. \n '...
            'Using default values.\n']);
    Data.maxeval=varargin{1};    %user defined maximum number of allowed function evaluations
    Data.surrogate_model=varargin{2}; %user defined surrogate model type
    Data.sampling_technique=varargin{3}; %user defined sampling technique
    Data.initial_design='LHS';  %if no initial experimental design strategy defined, use Latin hypercube
    Data.number_startpoints=2*(Data.dim+1); %default number of starting points
    Data.starting_point=[]; %by default no starting point defined
    Data.NumberNewSamples=1;  %by default only one point selected in every iteration
elseif length(varargin) < 5
    fprintf(['No number of starting points defined. \n '...
            'No starting points defined. \n '...
            'No number of new samples defined. \n '...
            'Using default value.\n ']);
    Data.maxeval=varargin{1};    %user defined maximum number of allowed function evaluations
    Data.surrogate_model=varargin{2}; %user defined surrogate model type
    Data.sampling_technique=varargin{3};  %user defined sampling technique
    Data.initial_design=varargin{4};     %user defined initial experimental design
    Data.number_startpoints=2*(Data.dim+1); %default number of starting points 
    Data.starting_point=[]; %by default no starting point defined
    Data.NumberNewSamples=1;  %by default only one point selected in every iteration
elseif length(varargin) < 6
    fprintf(['No starting points defined. \n '...
            'No number of new samples defined. \n '...
            'Using default value.\n ']);
    Data.maxeval=varargin{1};    %user defined maximum number of allowed function evaluations
    Data.surrogate_model=varargin{2}; %user defined surrogate model type
    Data.sampling_technique=varargin{3};  %user defined sampling technique
    Data.initial_design=varargin{4};    %user defined initial experimental design
    Data.number_startpoints=varargin{5};  %user defined number of starting points  
    Data.starting_point=[]; %by default no starting point defined
    Data.NumberNewSamples=1;  %by default only one point selected in every iteration
elseif length(varargin) < 7  %check if number of new samples in each iteration given
    fprintf(['No number of new samples defined. \n '...
            'Using default value.\n ']);
    Data.maxeval=varargin{1};     %user defined maximum number of allowed function evaluations
    Data.surrogate_model=varargin{2};  %user defined surrogate model type
    Data.sampling_technique=varargin{3};  %user defined sampling technique
    Data.initial_design=varargin{4};  %user defined initial experimental design
    Data.number_startpoints=varargin{5};   %user defined number of starting points  
    Data.starting_point=varargin{6}; %user defined starting point
    Data.NumberNewSamples=1;  %by dfault only one point selected in every iteration
else %all optional input arguments given
    Data.maxeval=varargin{1};    %user defined maximum number of allowed function evaluations
    Data.surrogate_model=varargin{2}; %user defined surrogate model type
    Data.sampling_technique=varargin{3}; %user defined sampling technique
    Data.initial_design=varargin{4}; %user defined initial experimental design
    Data.number_startpoints=varargin{5};   %user defined number of starting points  
    Data.starting_point=varargin{6};  %user defined starting point
    Data.NumberNewSamples=varargin{7}; %user defined number of points selected in every iteration
end
%check which input arguments were left blank
if isfield(Data,'number_startpoints') && isempty(Data.number_startpoints) %starting point(s)
    fprintf(['Using default number of points for initial experimental '...
        'design: (2x(dimension+1)=%d).\n'],2*(Data.dim+1))
    Data.number_startpoints=2*(Data.dim+1); %use default number of starting points 
end

%check if specified starting point is of correct dimension
if isfield(Data,'starting_point') && ~isempty(Data.starting_point)
    [m_sp,n_sp]=size(Data.starting_point);
    if n_sp ~= Data.dim
        if m_sp == Data.dim
            Data.starting_point=Data.starting_point';
        else
            error('Provided starting point has incorrect dimension.')
        end
    else
        if ~isempty(Data.integer) && any(abs(Data.starting_point(Data.integer)-round(Data.starting_point(Data.integer)))) >0
            error('One or more of the given starting points do not satisfy the integrality conditions')
        end
    end
end

%check if Number of new samples in each iteration is defined. 
if isfield(Data,'NumberNewSamples') &&  ~isempty(Data.NumberNewSamples)
    if Data.NumberNewSamples <= 0 ||  abs(Data.NumberNewSamples -...
            round(Data.NumberNewSamples)) >0
        error('NumberNewSamples  must be a positive integer.');
    end
elseif (isfield(Data,'NumberNewSamples') &&  isempty(Data.NumberNewSamples)) || ...
        ~isfield(Data,'NumberNewSamples') 
    fprintf(['No input for NumberNewSamples provided by user. \n '...
            'Using default value: NumberNewSamples=1. \n'])
    Data.NumberNewSamples=1;
end

%check if MATLAB Parallel Computing Toolbox is installed and start
%matlabpool
Data.parallel_eval=0;
v = ver; %check toolboxes installed
if ~any(strcmp('Parallel Computing Toolbox', {v.Name}))
    fprintf(['You do not have MATLAB Parallel Computing Toolbox. \n'...
        'Function evaluations will be done in serial. \n'])
else
    fprintf(['You do have MATLAB Parallel Computing Toolbox. \n '...
        'Evaluations will be done in parallel to reduce computation time.\n'])
    Data.parallel_eval=1;
    if matlabpool('size') == 0 %matlabpool not yet open
        matlabpool open
    end
end

%check if surrogate model type was left empty and assign default values if
%needed
if isempty(Data.surrogate_model) %surrogate model
    fprintf('Using default surrogate model: cubic RBF.\n')
    Data.surrogate_model='RBFcub';   %if no surrogate model defined, use cubic RBF
else % check if desired surrogate model is implemented; add own surrogate model to the list if needed
     if ~any(strcmp(Data.surrogate_model,{'RBFcub', 'RBFtps', 'RBFlin',...
             'POLYlin', 'POLYquad', 'POLYquadr', 'POLYcub', 'POLYcubr', ...
             'MARS', 'MIX_RcM','MIX_RcPc', 'MIX_RcPq','MIX_RcPcM','MIX_RcPqr','MIX_RcPcr'}))
         error('The surrogate model you want to use is not contained in the toolbox. Check spelling.')
     end
end
%check if sampling technique was left empty and assign default values if
%needed
if isempty(Data.sampling_technique) %sampling technique
    fprintf('Using default sampling strategy: CANDglob-global candidate point sampling.\n')
    Data.sampling_technique='CANDglob';  %if no sampling method definded, use candidate points
else %check if user defined sampling technique exists
    if ~any(strcmp(Data.sampling_technique,{'CANDloc','CANDglob','SurfMin'}))
        error('The sampling technique you want to use is not contained in the toolbox. Check spelling.')
    end  
end
%if user wants to use SurfMin sampling, check if optimization toolboxes are
%installed
if strcmp(Data.sampling_technique,'SurfMin') 
    v=ver;
    if ~isempty(Data.continuous) && ~any(strcmp('Optimization Toolbox', {v.Name}))
        disp(['You do not have MATLAB Optimization Toolbox. Cannot use option SurfMin. Will use option CANDglob.'])
        Data.sampling_technique='CANDglob'; 
    elseif ~isempty(Data.integer) && ~any(strcmp('Global Optimization Toolbox', {v.Name}))
        disp(['You do not have MATLAB Global Optimization Toolbox. Cannot use option SurfMin. Will use option CANDglob.'])
        Data.sampling_technique='CANDglob'; 
    end
end
%check if initial experimental design technique was left empty and assign 
%default values if needed               
if isempty(Data.initial_design) %initial experimental design
    fprintf('Using default initial experimental design strategy: LHS.\n')
    Data.initial_design='LHS';      %if no initial experimental design strategy defined, use Latin hypercube
else %check if given initial experimental design technique exists
    if ~any(strcmp(Data.initial_design,{'SLHD', 'LHS', 'CORNER'}))
        error('The experimental design strategy you want to use is not contained in the toolbox. Check spelling.')
    end
end
%check if maximum number of allowed evaluations was left empty and assign 
%default value if needed  
if isempty(Data.maxeval) %maximum number of function evaluations
    fprintf(['Using default number of maximal allowed function evaluations: 20xdimension=%d.\n'],20*Data.dim)
    Data.maxeval= 20*Data.dim;              
elseif ~isempty(Data.maxeval) && (Data.maxeval <= 0 || abs(Data.maxeval - round(Data.maxeval)) >0)
    error('The maximum number of allowed function evaluations must be a reasonably chosen positive integer.')
end

%check if there are enough starting points for using RBF models
if any(strcmp(Data.surrogate_model,{'RBFcub', 'RBFlin','RBFtps'})) 
    if isfield(Data,'starting_point') && (Data.number_startpoints+...
            size(Data.starting_point,1) <Data.dim+1)
        fprintf(['With RBF model I need at least Data.dim+1=%d starting points.\n'],Data.dim+1)
        Data.number_startpoints =Data.dim+1-size(Data.starting_point,1);
    elseif(Data.number_startpoints <Data.dim+1)  
        fprintf('With RBF model I need at least Data.dim+1=%d starting points.\n',Data.dim+1)
        Data.number_startpoints =Data.dim+1;
    end
end

%check if there are enough starting points for using linear polynomial regression 
%model
if any(strcmp(Data.surrogate_model,'POLYlin')) 
    if isfield(Data,'starting_point') && (Data.number_startpoints+...
            size(Data.starting_point,1) < 1+Data.dim)
        fprintf('For linear polynomial I need at least 1+Data.dim=%d starting points.\n',Data.dim+1)
        Data.number_startpoints =1+Data.dim-size(Data.starting_point,1);
    elseif (Data.number_startpoints <1+Data.dim)
        fprintf('For linear polynomial I need at least 1+Data.dim=%d starting points.\n',Data.dim+1)
        Data.number_startpoints =1+Data.dim;
    end
end

%check if there are enough starting points for using MARS model
if any(strcmp(Data.surrogate_model,'MARS')) 
    if isfield(Data,'starting_point') && (Data.number_startpoints+...
            size(Data.starting_point,1) < 1+Data.dim)
        fprintf('For the MARS model I need at least 1+Data.dim=%d starting points.\n',Data.dim+1)
        Data.number_startpoints =1+Data.dim-size(Data.starting_point,1);
    elseif (Data.number_startpoints <1+Data.dim)
        fprintf('For the MARS model I need at least 1+Data.dim=%d starting points.\n',Data.dim+1)
        Data.number_startpoints =1+Data.dim;
    end
end

%check if there are enough starting points for using full cubic regression 
%model
if any(strcmp(Data.surrogate_model,'POLYcub')) 
    if Data.dim > 2
        if isfield(Data,'starting_point') && (Data.number_startpoints+...
                size(Data.starting_point,1) < 1+3*Data.dim+ ...
                nchoosek(Data.dim,2)+ nchoosek(Data.dim,3))
            fprintf(['For full cubic polynomial I need at least\n 1+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)=%d starting points.\n',...
                'Otherwise the least squares problem is underdetermined.\n'],1+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3));
            Data.number_startpoints =1+3*Data.dim+ nchoosek(Data.dim,2)+ ...
                nchoosek(Data.dim,3)-size(Data.starting_point,1);
        elseif (Data.number_startpoints <1+3*Data.dim+ nchoosek(Data.dim,2)+ ...
                nchoosek(Data.dim,3)) 
            fprintf(['For full cubic polynomial I need at least\n 1+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)=%d starting points.\n',...
                'Otherwise the least squares problem is underdetermined.\n'],1+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3));
            Data.number_startpoints =1+3*Data.dim+ nchoosek(Data.dim,2)+nchoosek(Data.dim,3);
        end
        if strcmp(Data.initial_design,'CORNER') && Data.number_startpoints> 2^Data.dim+1
            error(['Bad combination of initial starting design and surrogate model.\n',...
                'POLYcub needs at least 1+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)=%d starting points,\n ',...
                'but CORNER allows only 2^Data.dim+1 =%d points.'],1+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3),2^Data.dim+1)
        end
    else
        error('Cannot use full cubic polynomial for problems with fewer than 3 (three) dimensions')
    end
end

%check if there are enough starting points for using reduced cubic 
%polynomial regression model
if any(strcmp(Data.surrogate_model,{'POLYcubr','MIX_RcPcr'})) 
    if isfield(Data,'starting_point') && (Data.number_startpoints+...
            size(Data.starting_point,1) < 1+3*Data.dim)
        fprintf(['For reduced cubic polynomial I need at least 1+3*Data.dim=%d starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'],1+3*Data.dim)
        Data.number_startpoints =1+3*Data.dim-size(Data.starting_point,1);
    elseif (Data.number_startpoints <1+3*Data.dim) 
        fprintf(['For reduced cubic polynomial I need at least 1+3*Data.dim=%d starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'],1+3*Data.dim);
        Data.number_startpoints =1+3*Data.dim;
    end
    if strcmp(Data.initial_design,'CORNER') && Data.number_startpoints> 2^Data.dim+1
        error(['Bad combination of initial starting design and surrogate model.\n',...
            'POLYcubr/MIX_RcPcr needs at least 1+3*Data.dim=%d starting points,\n ',...
            'but CORNER allows only 2^Data.dim+1 =%d points.'],1+3*Data.dim,2^Data.dim+1)
    end
end

%check if there are enough starting points for using full quadratic regression 
%model
if any(strcmp(Data.surrogate_model,'POLYquad')) 
    if Data.dim > 1
        if isfield(Data,'starting_point') && (Data.number_startpoints+...
                size(Data.starting_point,1) < 1+2*Data.dim+nchoosek(Data.dim,2))
            fprintf(['For full quadratic polynomial I need at least\n 1+2*Data.dim+ nchoosek(Data.dim,2)=%d starting points.\n',...
                'Otherwise the least squares problem is underdetermined.\n'],1+2*Data.dim+ nchoosek(Data.dim,2));
            Data.number_startpoints =1+2*Data.dim+ nchoosek(Data.dim,2) ...
               -size(Data.starting_point,1);
        elseif (Data.number_startpoints <1+2*Data.dim+ nchoosek(Data.dim,2)) 
            fprintf(['For full quadratic polynomial I need at least\n 1+2*Data.dim+ nchoosek(Data.dim,2)=%d starting points.\n',...
                'Otherwise the least squares problem is underdetermined.\n'],1+2*Data.dim+ nchoosek(Data.dim,2));
            Data.number_startpoints =1+2*Data.dim+ nchoosek(Data.dim,2);
        end
        if strcmp(Data.initial_design,'CORNER') && Data.number_startpoints> 2^Data.dim+1
            error(['Bad combination of initial starting design and surrogate model.\n',...
                'POLYquad needs at least 1+2*Data.dim+ nchoosek(Data.dim,2)=%d starting points,\n ',...
                'but CORNER allows only 2^Data.dim+1 =%d points.'],1+2*Data.dim+ nchoosek(Data.dim,2),2^Data.dim+1)
        end
    else
        error('Cannot use full quadratic polynomial for problems with fewer than 2 (two) dimensions')
    end
end
%check if there are enough starting points for using reduced quadratic 
%polynomial regression model
if any(strcmp(Data.surrogate_model,{'POLYquadr','MIX_RcPqr'})) 
    if isfield(Data,'starting_point') && (Data.number_startpoints+...
            size(Data.starting_point,1) < 1+2*Data.dim)
        fprintf(['For reduced quadratic polynomial I need at least 1+2*Data.dim=%d starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'],1+2*Data.dim)
        Data.number_startpoints =1+2*Data.dim-size(Data.starting_point,1);
    elseif (Data.number_startpoints <1+2*Data.dim) 
        fprintf(['For reduced quadratic polynomial I need at least 1+2*Data.dim=%d starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'],1+2*Data.dim)
        Data.number_startpoints =1+2*Data.dim;
    end
end
%check if there are enough starting points for using mixtures with full
%quadratic polynomial regression model
if any(strcmp(Data.surrogate_model,{'MIX_RcPq'})) 
    if Data.dim > 1
        if isfield(Data,'starting_point') && (Data.number_startpoints+...
                size(Data.starting_point,1) < 1+2*Data.dim+nchoosek(Data.dim,2)) %2+ because is mixture model
            fprintf(['For mixture with full quadratic polynomial I need at least\n 1+2*Data.dim+nchoosek(Data.dim,2)=%d starting points.\n',...
                'Otherwise the least squares problem is underdetermined.\n'],1+2*Data.dim+nchoosek(Data.dim,2))
            Data.number_startpoints =1+2*Data.dim+nchoosek(Data.dim,2)-size(Data.starting_point,1);
        elseif (Data.dim >=2)&& (Data.number_startpoints <1+2*Data.dim+nchoosek(Data.dim,2)) 
            fprintf(['For mixture with full quadratic polynomial I need at least\n 1+2*Data.dim+nchoosek(Data.dim,2)=%d starting points.\n',...
                'Otherwise the least squares problem is underdetermined.\n'],1+2*Data.dim+nchoosek(Data.dim,2))
            Data.number_startpoints =1+2*Data.dim+nchoosek(Data.dim,2);
        elseif Data.dim <2
            fprintf('For 1-dim problem there are no intersection terms of the form x1x2 possible.\n');
        end
        if strcmp(Data.initial_design,'CORNER') && Data.number_startpoints> 2^Data.dim+1
            error(['Bad combination of initial starting design and surrogate model.\n',...
            'MIX_RcPq needs at least 1+2*Data.dim+nchoosek(Data.dim,2)=%d starting points,\n ',...
            'but CORNER allows only 2^Data.dim+1 =%d points.'],1+2*Data.dim+nchoosek(Data.dim,2),2^Data.dim+1)
        end
    else
        error('Cannot use mixture with full quadratic polynomial for problems with fewer than 2 (two) dimensions')
    end
end

%check if there are enough starting points for using mixtures with full
%cubic polynomial regression model
if any(strcmp(Data.surrogate_model,{'MIX_RcPc', 'MIX_RcPcM'})) 
    if Data.dim > 2
        if isfield(Data,'starting_point') && (Data.number_startpoints+...
                size(Data.starting_point,1) < 2+3*Data.dim+ ...
                nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)) %2+ because is mixture model
            fprintf(['For mixture with full cubic polynomial I need at least\n 2+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)=%d starting points.\n',...
                'Otherwise the least squares problem is underdetermined.\n'],2+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3))
            Data.number_startpoints =2+3*Data.dim+ nchoosek(Data.dim,2)+...
                nchoosek(Data.dim,3)-size(Data.starting_point,1);
        elseif (Data.dim >=3)&& (Data.number_startpoints <2+3*Data.dim+ ...
                nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)) 
            fprintf(['For mixture with full cubic polynomial I need at least\n 2+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)=%d starting points.\n',...
                'Otherwise the least squares problem is underdetermined.\n'],2+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3))
            Data.number_startpoints =2+3*Data.dim+ nchoosek(Data.dim,2)+...
                nchoosek(Data.dim,3);
        elseif Data.dim <3
            fprintf('For 2-dim problem there are no intersection terms of the form x1x2x3 possible.\n');
        end
        if strcmp(Data.initial_design,'CORNER') && Data.number_startpoints> 2^Data.dim+1
            error(['Bad combination of initial starting design and surrogate model.\n',...
            'MIX_RcPc/MIX_RcPcM needs at least 2+3*Data.dim+nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)=%d starting points,\n ',...
            'but CORNER allows only 2^Data.dim+1 =%d points.'],2+3*Data.dim+ ...
                nchoosek(Data.dim,2)+ nchoosek(Data.dim,3),2^Data.dim+1)
        end
    else
        error('Cannot use mixture with full cubic polynomial for problems with fewer than 3 (three) dimensions')
    end
end

%check if there are enough starting points for using mixtures with MARS or
%radial basis functions
if any(strcmp(Data.surrogate_model,{'MIX_RcM'})) 
    if isfield(Data,'starting_point') && (Data.number_startpoints+...
            size(Data.starting_point,1) < Data.dim+2) %2+ because is mixture model
        fprintf(['For the chosen mixture I need at least 2+Data.dim=%d starting points.\n'],2+Data.dim)
        Data.number_startpoints =2+Data.dim-size(Data.starting_point,1);
    elseif (Data.number_startpoints <2+Data.dim) 
        fprintf(['For the chosen mixture I need at least 2+Data.dim=%d starting points.\n'],2+Data.dim)
        Data.number_startpoints =2+Data.dim;
    end
end

%check how many integer variables are 0-1 variables
if ~isempty(Data.integer)
    n01_integer = length(find(Data.xup(Data.integer)-Data.xlow(Data.integer)==1));
    if ((strcmp(Data.surrogate_model(1:4),'POLY')) ||  any(strcmp(Data.surrogate_model,{'MIX_RcPc','MIX_RcPq','MIX_RcPcM'}))) &&...
            n01_integer  > 0
      error(['Using a polynomial regression model or ensemble containing a '...
          'polynomial regression model is not possible when there are '...
          'integers that are binary. Try an RBF model.'])
    end
end

%check if there are enough starting points for using SLHD design or CORNER
if strcmp(Data.initial_design,'SLHD') && exist('starting_point','var')
    if (Data.number_startpoints+size(Data.starting_point,1) <2*Data.dim)
        fprintf(['With SLHD I need to use at least 2*Data.dim=%d starting points.\n'],2*Data.dim)
        Data.number_startpoints = 2*Data.dim-size(Data.starting_point,1);
    end
elseif strcmp(Data.initial_design,'SLHD') 
    if (Data.number_startpoints < 2*Data.dim)
        fprintf(['With SLHD I need to use at least 2*Data.dim=%d starting points.\n'],2*Data.dim)
        Data.number_startpoints = 2*Data.dim;
    end
elseif strcmp(Data.initial_design,'CORNER') && exist('starting_point','var')
    if (Data.number_startpoints+size(Data.starting_point,1) > 2^Data.dim+1)
        fprintf(['With CORNER as initial design strategy, I can have at most 1+2^Data.dim=%d points.\n',...
            'Using the first %d points to create the initial experimental design.\n'],1+2^Data.dim,1+2^Data.dim)
        Data.number_startpoints = min(2^Data.dim+1,Data.number_startpoints)-size(Data.starting_point,1);
    end
elseif strcmp(Data.initial_design,'CORNER') 
    if (Data.number_startpoints > 2^Data.dim+1)
        fprintf(['With CORNER as initial design strategy, I can have at most 1+2^Data.dim=%d points.\n',...
            'Using only %d points to create the initial experimental design.\n'],1+2^Data.dim,1+2^Data.dim)
        Data.number_startpoints = 2^Data.dim+1;
    end    
end

%check if the number of defined or necessary starting points exceeds budget
%of allowable function evaluations
if Data.number_startpoints > Data.maxeval
    error(['The number of starting points required to fit the surrogate '...
        'model exceeds the total number of allowed function evaluations. '...
        'Try a different response surface. '])
end
%% -----END INPUT CHECK-------
end %function
