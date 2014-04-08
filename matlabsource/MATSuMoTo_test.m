function MATSuMoTo_test

disp('Testing continuous problem')
MATSuMoTo('datainput_hartman3');
disp('Continuous problem finished successfully')

disp('Testing pure integer problem')
MATSuMoTo('datainput_convex_I');
disp('Pure integer problem finished successfully')

disp('Testing pure integer problem')
MATSuMoTo('datainput_nvs09_MI',60);
disp('Mixed integer problem finished successfully')

MATSuMoTo('datainput_Powell24',500,[],[],'CORNER');

data_file = 'datainput_hartman3'; %name of data file
maxeval = 300; %maximum number of allowed function evaluations
surogate_model = 'RBFtps'; %selected surrogate model type
sampling_technique = 'CANDglob'; %global randomized sampling strategy
initial_design = 'LHS'; %Matlab's lhsdesign.m as initial design
number_startpoints = 15; %15 points in the initial experimental design
starting_point = []; %no user-specified points to be added to the initial design
NumberNewSamples = 1; %1 new point is selected in each iteration
MATSuMoTo(data_file,maxeval,surogate_model,sampling_technique,...
    initial_design,number_startpoints,starting_point,NumberNewSamples);


data_file = 'datainput_Shekel7'; %name of data file
maxeval = 500; %maximum number of allowed function evaluations
surogate_model = 'MIX_RcPcr'; %selected response surface
sampling_technique = 'CANDglob'; %global randomized sampling strategy
initial_design = 'SLHD'; %symmetric Latin hypercube design as initial design
number_startpoints = []; %default number of points for initial experimental design
starting_point = []; %no user-specified points to be added to the initial design
NumberNewSamples = 1; %1 new point is selected in each iteration
MATSuMoTo(data_file,maxeval,surogate_model,sampling_technique,...
    initial_design,number_startpoints,starting_point,NumberNewSamples);