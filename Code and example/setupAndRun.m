function [Results,DatStore,Misc] = setupAndRun(Misc,passiveParm,funWeights)
%% Input data
% get time from IK solution
for i=1:length(Misc.IKfile)
    IK = importdata(Misc.IKfile{i});
    time(i,:) = [IK.data(1,1) IK.data(end,1)];
end

side = Misc.side;
Misc.DofNames_Input={['ankle_angle_' side],['knee_angle_' side],['hip_flexion_' side],['hip_adduction_' side],['hip_rotation_' side]};  
Misc.kT = []; % default 35      
%% Settings
% % select muscles
modlist={'bflh'  'semiten'  'vaslat'  'vasmed'  'recfem'  'tfl'  'glmed3'  'glmax1'  'tibant'  'gaslat'  'gasmed'  'soleus'  'vasint'};

% Provide the correct headers in case you EMG file has not the same
% headers as the muscle names in OpenSim (leave empty when you don't want
% to use this)
Misc.EMGheaders = {'time' [modlist{1} '_' side] [modlist{2} '_' side] [modlist{3} '_' side] [modlist{4} '_' side] [modlist{5} '_' side] [modlist{6} '_' side]...
                          [modlist{7} '_' side] [modlist{8} '_' side] [modlist{9} '_' side] [modlist{10} '_' side] [modlist{11} '_' side] [modlist{12} '_' side] [modlist{13} '_' side]};

% channels you want to use for EMG constraints
Misc.EMGSelection = {[     modlist{1} '_' side] [modlist{2} '_' side] [modlist{3} '_' side] [modlist{4} '_' side] ...
                                                                      [modlist{9} '_' side] [modlist{10} '_' side] [modlist{11} '_' side] [modlist{12} '_' side] [modlist{13} '_' side]};
                                              
% Settings for estimating tendon stiffness
if strcmp(Misc.SimulationLabel,'emg_ten')
    Misc.Estimate_TendonStiffness = {[modlist{3} '_' side] [modlist{4} '_' side] [modlist{13} '_' side]    [modlist{10} '_' side] [modlist{11} '_' side] [modlist{12} '_' side]}; % Names of muscles of which tendon stifness is estimated
    Misc.lb_kT_scaling = 0.1; % Lower bound for scaling generic tendon stiffness
    Misc.ub_kT_scaling = 1.2; % Upper bound for scaling generic tendon stiffness
    Misc.Coupled_TendonStiffness  = {[modlist{3} '_' side] [modlist{4} '_' side] [modlist{13} '_' side];    [modlist{10} '_' side] [modlist{11} '_' side] [modlist{12} '_' side]}; % Couple muscles that should have equal tendon stifness -> plantarflexors and vasti
else 
    % choose the nonlinear tendon stiffness. Use the calibrated parameters from emg_ten to reproduce effort_ten
    PF_kT=35; VA_kT=35; % Choose any value for the normalized tendon stiffness
    Misc.Set_kT_ByName = {['soleus_' side],PF_kT;  ['gasmed_' side],PF_kT; ['gaslat_' side],PF_kT;    ['vasint_' side],VA_kT;  ['vaslat_' side], VA_kT;   ['vasmed_' side],VA_kT};
end

Misc.Coupled_EMG_scaling={[modlist{3} '_' side] [modlist{4} '_' side]; [modlist{3} '_' side] [modlist{13} '_' side]; [modlist{10} '_' side] [modlist{11} '_' side]};

Misc.EMGbounds      = [-0.01 0.01]; % upper and lower bound for difference between simulated and measured muscle activity
Misc.BoundsScaleEMG = [0.05 2.50];  % maximal value to scale EMG

% Assign passive parameters
Misc.PassiveParam_kpe= passiveParm(1,:);
Misc.PassiveParam_ksf= passiveParm(2,:);
Misc.PassiveParam_ke0= passiveParm(3,:);

% name output
Misc.OutName = ['simulation_']; 
Misc.MuscleNames_Input = []; % all muscle in the model

% run muscle analysis
Misc.RunAnalysis = true;

% weights
Misc.wEMG  = funWeights(1);
Misc.wAct  = funWeights(2);
Misc.wTres = funWeights(3);
Misc.wVm   = funWeights(4);
%% Run muscle tendon estimator:
[Results,DatStore,Misc] = solveMuscleRedundancy(time,Misc);
end