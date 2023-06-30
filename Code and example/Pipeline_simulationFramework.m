%% PIPELINE
% 
% This script serves as an example to run the pipeline. 
% We examplify three of the simulation workflows in our study: effort_gen, 
% effort_cal, and emg_ten. Simulation workflow "effort_ten" is easily reproduced 
% by changing the tendon stiffness in setAndRun function.
%% DESCRIPTION
%
% You can set up the type of simulation framework you want to run
% effort_gen -> Minimal muscle effort with generic passive force
% effort_cal -> Minimal muscle effort with calibrated passive force
% emg_ten    -> EMG-informed with calibrated passive force 
%
% Each simulation workflow modify the following features in the muscle
% redundacy solver (MRS)
%     Misc.MRSBool       -> MRS Bool: Select if you want to run the generic muscle redundancy solver
%     Misc.EMGconstr     -> Constraint Bool: Select if you want to constraint with EMGs
%     Misc.BoolParamOpt  -> Parameter optimization Bool: Select if you want to optimize parameters
%     Misc.PlotBool      -> Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
%     Misc.ValidationBool-> Validation Bool: Select if you want to run the
%     muscleredundacy solver with the optimized parameters from the
%     EMG-constraint simulation (we do not use this feature as it is not
%     required)
%     Misc.PassForceMode  -> Choose to use generic [genParam] or calibrated [calParam] passive parameters
%     funWeights -> Set the weights for the objective function: [EMG_weight, Act_weight, Res_Weight, Vm_weight].
%
%% EXPERIMENTAL DATA
clear all;
MainPath = pwd;
DataPath = fullfile(MainPath,'Data');

Misc.model_path = fullfile(DataPath,'1_scaledHMK_rajagopal2022_MIF2.osim');
Misc.IKfile  = strcat([DataPath '\'],[{'IK_sub1_slow.mot'} {'IK_sub1_normal.mot'} {'IK_sub1_fast.mot'}]);
Misc.IDfile  = strcat([DataPath '\'],[{'ID_sub1_slow.sto'} {'ID_sub1_normal.sto'} {'ID_sub1_fast.sto'}]);
Misc.EMGfile = strcat([DataPath '\'],[{'EMG_sub1_slow.mot'} {'EMG_sub1_normal.mot'} {'EMG_sub1_fast.mot'}]);
Misc.side    = 'r';
Misc.PassForcePath = fullfile(MainPath,'CalibrationOutcome','passiveMoments.mat'); % add the path where the passive parameters are located

%% SELECT SIMULATION WORKFLOW
framework_effort_gen = 1;
framework_effort_cal = 0;
framework_emg_ten    = 0;

if framework_effort_gen==1 % EFFORT-GEN 
    Misc.SimulationLabel= 'effort_gen';
    Misc.MRSBool       = 1;
    Misc.EMGconstr     = 0;
    Misc.BoolParamOpt  = 0;
    Misc.PlotBool      = 1;
    Misc.ValidationBool= 0;
    Misc.PassForceMode  = {'genParam'};
    funWeights = [0,1,1000,0.01];
end
if framework_effort_cal==1 % EFFORT-CAL
    Misc.SimulationLabel= 'effort_cal';
    Misc.MRSBool       = 1;
    Misc.EMGconstr     = 0;
    Misc.BoolParamOpt  = 0;
    Misc.PlotBool      = 1;
    Misc.ValidationBool= 0;
    Misc.PassForceMode  = {'calParam'};
    funWeights = [0,1,1000,0.01];
end
if framework_emg_ten==1 % EMG-TEN
    Misc.SimulationLabel= 'emg_ten';
    Misc.MRSBool       = 0;
    Misc.EMGconstr     = 1;
    Misc.BoolParamOpt  = 1;
    Misc.PlotBool      = 1;
    Misc.ValidationBool= 0;
    Misc.PassForceMode = {'calParam'};
    funWeights = [10,0.1,1000,0.01];
end
%% SETUP AND RUN SIMULATION
Misc.OutPath = fullfile(MainPath,'Results',Misc.SimulationLabel);
[passiveParm]= loadPassiveParam(Misc.PassForceMode,Misc.PassForcePath);
[Results,DatStore,Misc]= setupAndRun(Misc,passiveParm,funWeights);