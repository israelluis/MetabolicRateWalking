%% PASSIVE MUSCLE FORCES
%
% This script serves to run the calibration of passive forces. 
%% DESCRIPTION
% We first setup the information about the model, and the kinematics and
% dynamics that describe the relationship between passive torque and angle
% at the ankle, knee and hip. Then we compute the muscle moment arms,
% muscle-tendon lengths, and initial guess for the optimization,
% formulate and run the optimization routine. Later, the store the
% results in a mat file. Finally, the results are shown: Table and graph. 
%
%% Add paths and casadi
clear all; clc; close all;
MainPath = pwd;
DataPath = fullfile(MainPath,'Data');

addpath(genpath(MainPath));
addpath(genpath('C:\Users\israe\Dropbox\PhDWork\Simulation\code\MuscleRedundancySolver-master\Casadi_installed')); % <- Add your folder where casadi is installed

import casadi.* %x = MX.sym('x')
import org.opensim.modeling.*
%% INPUT INFORMATION

model_conf={'rajagopal'}; % 2392 (knee bending is negative) / rajagopal
side_sel  ='r';           % choose the leg side

Misc.model_path  = fullfile(DataPath,'Rajagopal2015_passiveCal_hipAbdMoved.osim');
Misc.OutPath     = fullfile(MainPath,'CalibrationOutcome');                         % folder to store results

% data based on Silder2007 (IDENTIFICATION OF PASSIVE ELASTIC JOINT MOMENT-ANGLE RELATIONSHIPS IN THE LOWER EXTREMITY, )
ank_file_name=fullfile(DataPath,['TrialAnk_' side_sel '_kne0hip0_kneBen15hip0_kneBen60hip0_' model_conf{:} 'Model']);
hip_file_name=fullfile(DataPath,['TrialHip_' side_sel '_kneBen15ank0_kneBen60ank0_' model_conf{:} 'Model']);
kne_file_name=fullfile(DataPath,['TrialKne_' side_sel '_ankDor20hip0_ankPla15hip0_ankDor20hipExt15_' model_conf{:} 'Model']);

Misc.IKfile = {[ank_file_name '.mot'], [hip_file_name '.mot'], [kne_file_name '.mot']};
Misc.IDfile = {[ank_file_name '.sto'], [hip_file_name '.sto'], [kne_file_name '.sto']};
           
% name output
Misc.OutName = 'ParmEst';

% input dofs
Misc.DofNames_Input={{['ankle_angle_' side_sel]};{['hip_flexion_' side_sel]};{['knee_angle_' side_sel]}};    % select the DOFs you want to include in the optimization

% adapt the stiffness of the achilles tendon (optional input argument)
% 2392
if strcmp(model_conf{:},{'2392'})
    Misc.Set_kT_ByName = {['soleus_' side_sel],20;  ['med_gas_' side_sel],20;    ['lat_gas_' side_sel],20};
end

% rajagopal
if strcmp(model_conf{:},{'rajagopal'})
    Misc.Set_kT_ByName = {['soleus_' side_sel],20;  ['gasmed_' side_sel],20;    ['gaslat_' side_sel],20};
end

time = [0 10; 0 10; 0 10];

% run muscle analysis
Misc.RunAnalysis = true;

%% Run muscle tendon estimator:

% No filter of data (as we assume that velocity = 0 )
Misc.f_order_dM = NaN;
Misc.f_order_lMT= NaN;
Misc.f_order_IK = NaN;
Misc.f_order_ID = NaN;
Misc.use_filter = 0;

% update default settings
Misc = DefaultSettings(Misc);

% read the muscle properties
[Misc] = getMuscleProperties(Misc.model_path,Misc);

% number of motion trials
Misc.nTrials = length(Misc.IKfile);

%check if we have to adapt the start and end time so that it corresponds to
%time frames in the IK solution
[time] = Check_TimeIndices(Misc,time);
Misc.time=time;

%% Perform muscle analysis for all trials
DatStore = struct;
MuscleAnalysisPath=fullfile(Misc.OutPath,'MuscleAnalysis');
Misc.MuscleAnalysisPath=MuscleAnalysisPath;
if ~exist(MuscleAnalysisPath,'dir')
    mkdir(MuscleAnalysisPath);
end
for i = 1:Misc.nTrials
    % select the IK and ID file
    IK_path_trial = Misc.IKfile{i};
    % Run muscle analysis
    if Misc.RunAnalysis
        disp('MuscleAnalysis Running .....');
        OpenSim_Muscle_Analysis(IK_path_trial,Misc.model_path,MuscleAnalysisPath,[time(i,1) time(i,end)],Misc.DofNames_Input{i,1})
        disp('MuscleAnalysis Finished');
    end
end

%% Extract muscle information
% Get number of degrees of freedom (dofs), muscle-tendon lengths, moment
% arms, stiffness and shift for the selected muscles.
for trial = 1:Misc.nTrials
    [~,Misc.MAtrialName{trial},~]=fileparts(Misc.IKfile{trial});
end

% select muscles with a moment arm for the input dofs
Misc = getMuscles4DOFS(Misc);

% get IK, ID, muscle analysis data
[Misc,DatStore] = getMuscleInfo(Misc,DatStore);

% display warnings in muscle selection
[Misc] = Warnings_MuscleNames(DatStore,Misc);

% get indexes of the muscles for which optimal fiber length, tendon stiffness are estimated
[DatStore,Misc] = GetIndices(DatStore,Misc);

% get the EMG and ultrasound information
[Misc,DatStore] = GetEMGInfo(Misc,DatStore);
[Misc,DatStore] = GetUSInfo(Misc,DatStore);

% get the number of muscles in a vector
NMuscles = zeros(Misc.nTrials,1);
for trial = 1:Misc.nTrials
    NMuscles(trial) = DatStore(trial).nMuscles;
end

%% Set solver
% I'm using ipopt, but there are probably better solvers for this simple
% problem
% Create an NLP solver
SolverSetup.nlp.solver = 'ipopt';
SolverSetup.derivatives.derivativelevel = 'second';
SolverSetup.optionssol.ipopt.nlp_scaling_method = 'gradient-based';
SolverSetup.optionssol.ipopt.linear_solver = 'mumps';
SolverSetup.optionssol.ipopt.tol = 1e-6;
SolverSetup.optionssol.ipopt.max_iter = 10000;

%% Initial guess (IG)
% initial guess is based on rigid tendon (and ignores pennation angle)
plot_initial_fiber_IG=0; % Boolean: Select if you want to see the initial guess of the fiber lengths

for trial = 1:Misc.nTrials
    nMuscles = length(DatStore(trial).MuscleNames);
    Nfr = length(DatStore(trial).time);
    lTs = Misc.lTs(Misc.idx_allMuscleList{trial})';
    IG(trial).lM_projected = zeros(nMuscles,Nfr);
    IG(trial).lMtilde = zeros(nMuscles,Nfr);
    for k=1:Nfr
        lMT = DatStore(trial).LMT(k,:)';
        lMGuess = lMT-lTs;
        lMo = Misc.params(2,Misc.idx_allMuscleList{trial})';
        alphao = Misc.params(4,Misc.idx_allMuscleList{trial})';
        w = lMo.*sin(alphao);
        IG(trial).lM_projected(:,k) = sqrt((lMGuess.^2 - w.^2));
        IG(trial).lMtilde(:,k) = lMGuess./lMo;
    end
    
    if plot_initial_fiber_IG==1
    figure();
    plot(IG(trial).lMtilde');
    xlabel('frames');
    ylabel('rigid tendon estimate fiber length');
    end
end

%% FORMULATE THE OPTIMIZATION PROBLEM
% we mimize difference between inverse dynamic moments and moments
% generated by muscles. We account for tendon compliance here and solve
% for equilibrium between tendon force and the force in the passive element
% in the muscle fiber

import casadi.*
opti    = casadi.Opti();    % create opti structure
J = 0;
for trial=1:Misc.nTrials
    % get number of muscles and frames
    nMuscles = length(DatStore(trial).MuscleNames);
    Nfr = length(DatStore(trial).time);

    % optimization variables
    TR(trial).kpe = opti.variable(nMuscles,1);
    TR(trial).ksf = opti.variable(nMuscles,1);

    % bounds
    ke0=0.6.*ones(nMuscles,1);
    opti.subject_to(3.0 < TR(trial).kpe < 5.0); % 4 by default
    opti.subject_to(0.8 < TR(trial).ksf < 1.2); % 1 by default
    
    % set initial gues
    opti.set_initial(TR(trial).kpe, 4);
    opti.set_initial(TR(trial).ksf, 1);

    % projected fiber length as optimization variable
    TR(trial).lM_projected = opti.variable(nMuscles,Nfr);
    opti.subject_to(1e-4 < TR(trial).lM_projected(:)); % We impose that projected muscle fiber length has strict positive length
    opti.set_initial(TR(trial).lM_projected,IG(trial).lM_projected);

    % create reserve actuators
    TR(trial).aTk = opti.variable(DatStore(trial).nDOF,Nfr);
    opti.subject_to(-1< TR(trial).aTk < 1)

    % lMtilde as optimization variable
    TR(trial).lMtilde = opti.variable(nMuscles,Nfr);
    opti.subject_to(0.05 < TR(trial).lMtilde < 2);
    opti.set_initial(TR(trial).lMtilde,IG(trial).lMtilde);

    % select muscle properties
    MuscProperties.params = Misc.params(:,Misc.idx_allMuscleList{trial})';
    MuscProperties.kT = Misc.kT(:,Misc.idx_allMuscleList{trial}')';
    MuscProperties.shift = Misc.shift(:,Misc.idx_allMuscleList{trial}')';

%     params_casadi =MX(nMuscles,8);
    params_casadi =[MuscProperties.params(:,1:5) TR(trial).kpe TR(trial).ksf ke0];
    
    % constraint on projected fiber length
    lMo = Misc.params(2,Misc.idx_allMuscleList{trial})';
    alphao = Misc.params(4,Misc.idx_allMuscleList{trial})';
    lM = TR(trial).lMtilde.*lMo;
    w = lMo.*sin(alphao);
    opti.subject_to(lM.^2 - w.^2 == TR(trial).lM_projected.^2);

    % find passive
    TR(trial).TsimVect = MX(Nfr,DatStore(trial).nDOF);
    for k=1:Nfr
        % solve muscle equilibrium
        [err, FTk] = ForceEquilibrium_lMtildeState_optPassive(0,TR(trial).lMtilde(:,k),...
            0,TR(trial).lM_projected(:,k),DatStore(trial).LMT(k,:)',params_casadi,...
            MuscProperties.kT,MuscProperties.shift);

        % Moment constraints
        for dof = 1:DatStore(trial).nDOF
            T_exp = DatStore(trial).T_exp(k,dof); % + 1 due to the time col
            T_sim = FTk'*squeeze(DatStore(trial).dM(k,dof,:)) + Misc.Topt*TR(trial).aTk(dof,k);
            opti.subject_to(T_exp - T_sim == 0);
            TR(trial).TsimVect(k,dof) = FTk'*squeeze(DatStore(trial).dM(k,dof,:));
        end
        % impose equilibrium
        opti.subject_to(err == 0);
    end
    J = J + sumsqr(TR(trial).aTk);
end

opti.minimize(J);
opti.solver(SolverSetup.nlp.solver,SolverSetup.optionssol);
optiDef = opti.copy(); % copy for solution with default params

% impose contraints such that bi-articular muscles that appear in multiple
% phases in the optimization problem have always the same parameters
Trial1Copies = nan(10,3);
ct = 1;
for iM = 1:length(DatStore(1).MuscleNames)
    iTrial2 = find(strcmp(DatStore(2).MuscleNames,DatStore(1).MuscleNames{iM}));
    iTrial3 = find(strcmp(DatStore(3).MuscleNames,DatStore(1).MuscleNames{iM}));
    if ~isempty(iTrial2)
        Trial1Copies(ct,1) = iM;
        Trial1Copies(ct,2) = iTrial2;
        Trial1Copies(ct,3) = 2;
        ct = ct+1;
    end
    if ~isempty(iTrial3)
        Trial1Copies(ct,1) = iM;
        Trial1Copies(ct,2) = iTrial3;
        Trial1Copies(ct,3) = 3;
        ct = ct+1;
    end
end
Trial1Copies(ct:end,:) = [];

Trial2Copies = nan(10,3);
ct = 1;
for iM = 1:length(DatStore(2).MuscleNames)
    iTrial3 = find(strcmp(DatStore(3).MuscleNames,DatStore(2).MuscleNames{iM}));   
    if ~isempty(iTrial3)
        Trial2Copies(ct,1) = iM;
        Trial2Copies(ct,2) = iTrial3;
        Trial2Copies(ct,3) = 3;
        ct = ct+1;
    end
end
Trial2Copies(ct:end,:) = [];

% set equality constraints
for i=1:length(Trial1Copies(:,1))
    % equality constraint
    iM_tr1 = Trial1Copies(i,1);
    iM_trx = Trial1Copies(i,2);
    trx = Trial1Copies(i,3);
    opti.subject_to(TR(1).kpe(iM_tr1,1) == TR(trx).kpe(iM_trx));
    opti.subject_to(TR(1).ksf(iM_tr1,1) == TR(trx).ksf(iM_trx));
end
for i=1:length(Trial2Copies(:,1))
    % equality constraint
    iM_tr2 = Trial2Copies(i,1);
    iM_trx = Trial2Copies(i,2);
    trx = Trial2Copies(i,3);
    opti.subject_to(TR(2).kpe(iM_tr2,1) == TR(trx).kpe(iM_trx));
    opti.subject_to(TR(2).ksf(iM_tr2,1) == TR(trx).ksf(iM_trx));
end

%% Additional constraint 
% muscle groups have the same values (since we for example cannot descriminate between parameters of the vastus
% lateral, medialis and intermdius we impose that the parameters should be the same)
% muscle_coupling_passive=cell(1,Misc.nTrials);
for trial=1:Misc.nTrials
    for i=1:length(Misc.DofNames_Input{trial})
        dof_name=Misc.DofNames_Input{trial}{i}(1:end-2);  %degree of freedom without side
        side    =Misc.DofNames_Input{trial}{i}(end);%side only
        [muscle_coupling_passive{trial}]= muscleCouplingPassive(model_conf,dof_name,side);
    end
end

for trial=1:Misc.nTrials
    [row,col]=size(muscle_coupling_passive{trial});
    MuscleNames=DatStore(trial).MuscleNames;
    for i=1:row
        index1=find(contains(MuscleNames,muscle_coupling_passive{trial}(i,1)));
        index2=find(contains(MuscleNames,muscle_coupling_passive{trial}(i,2)));
        opti.subject_to(TR(trial).kpe(index1) == TR(trial).kpe(index2));
        opti.subject_to(TR(trial).ksf(index1) == TR(trial).ksf(index2));
        disp(['Constraint on muscles ' DatStore(trial).MuscleNames(index1) DatStore(trial).MuscleNames(index2)])
    end
end
%% CALIBRATION COMPUTATION
% solve the optimization problem
Sol = opti.solve();
%% GENERIC COMPUTATION
import casadi.*

% solve with default values
for trial=1:Misc.nTrials
    % equality constraints for parameters of passive force length
    optiDef.subject_to( TR(trial).kpe == 4);
    optiDef.subject_to( TR(trial).ksf == 1 );
    nMuscles = length(DatStore(trial).MuscleNames);
    Nfr = length(DatStore(trial).time);
    % Loop over frames just to compute the torque generqted by the muscles
    % (this is actually postprocessing to create the matrix TsimVect)
    MuscProperties.params = Misc.params(:,Misc.idx_allMuscleList{trial})';
    MuscProperties.kT = Misc.kT(:,Misc.idx_allMuscleList{trial}')';
    MuscProperties.shift = Misc.shift(:,Misc.idx_allMuscleList{trial}')';
    
    ke0=0.6.*ones(nMuscles,1);
    params_casadi =[MuscProperties.params(:,1:5) TR(trial).kpe TR(trial).ksf ke0];
        
    TRDef(trial).TsimVect = MX(Nfr,DatStore(trial).nDOF);
    for k=1:Nfr
        [err, FTk] = ForceEquilibrium_lMtildeState_optPassive(0,TR(trial).lMtilde(:,k),...
            0,TR(trial).lM_projected(:,k),DatStore(trial).LMT(k,:)',params_casadi,...
            MuscProperties.kT,MuscProperties.shift);
       
        for dof = 1:DatStore(trial).nDOF
            T_sim = FTk'*squeeze(DatStore(trial).dM(k,dof,:)) + Misc.Topt*TR(trial).aTk(dof,k);
            TRDef(trial).TsimVect(k,dof) = FTk'*squeeze(DatStore(trial).dM(k,dof,:));
        end
    end
end
SolDefault = optiDef.solve();
%% Get solutions
calParam=struct.empty;
genParam=struct.empty;

for trial=1:Misc.nTrials
    % solution with parameter optimization
    calParam(trial).kpe = Sol.value(TR(trial).kpe);
    calParam(trial).ksf = Sol.value(TR(trial).ksf);
    calParam(trial).aTk = Sol.value(TR(trial).aTk).*Misc.Topt;
    calParam(trial).Tsim = Sol.value(TR(trial).TsimVect);
    % solution with default parameters
    genParam(trial).kpe = SolDefault.value(TR(trial).kpe);
    genParam(trial).ksf = SolDefault.value(TR(trial).ksf);
    genParam(trial).aTk = SolDefault.value(TR(trial).aTk).*Misc.Topt;
    genParam(trial).Tsim = SolDefault.value(TRDef(trial).TsimVect);
end

%% Save results
import org.opensim.modeling.*

model      = Model(Misc.model_path);
NumMuscle  = model.getMuscles().getSize;
MuscleNames=cell(1,NumMuscle);

for i=1:NumMuscle
    MuscleNames{i} = char(model.getMuscles().get(i-1).getName());
end

muscle_kpe_Cal=NaN(1,NumMuscle);    muscle_ksf_Cal=NaN(1,NumMuscle);
muscle_kpe_Def=NaN(1,NumMuscle);    muscle_ksf_Def=NaN(1,NumMuscle);

for trial=1:Misc.nTrials
    names=DatStore(trial).MuscleNames;
    indexes=find(contains(MuscleNames,names));
    muscle_kpe_Cal(indexes)=calParam(trial).kpe;
    muscle_ksf_Cal(indexes)=calParam(trial).ksf;
    muscle_kpe_Def(indexes)=genParam(trial).kpe;
    muscle_ksf_Def(indexes)=genParam(trial).ksf;    
end

Results.calParam=calParam;
Results.genParam=genParam;

Results.onlyParam.calParam.kpe=muscle_kpe_Cal;
Results.onlyParam.calParam.ksf=muscle_ksf_Cal;
Results.onlyParam.genParam.kpe=muscle_kpe_Def;
Results.onlyParam.genParam.ksf=muscle_ksf_Def;
Results.onlyParam.MuscleNames=MuscleNames;

sim.DatStore=DatStore;
sim.Misc=Misc;
sim.Results=Results;

save(fullfile(Misc.OutPath,'passiveMoments.mat'),'sim')
%% Plot variables as a list
kpe_string=cellstr(num2str(muscle_kpe_Cal(1:40)','%4.1f'));
ksf_string=cellstr(num2str(muscle_ksf_Cal(1:40)','%4.1f'));
disp( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' );
disp( '                Summary of outcomes' );
disp( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' );
disp( [[{'muscle name'} {'kpe value'} {'s0 value'}]; [MuscleNames(1:40)' kpe_string ksf_string]]);

%% Plot outcome as a time-series
fig=figure(1); clf; set(gcf,'color','w','Visible','on');

trial_sel =[3 2 3]; % # of trial per file
trial_list=[1 3 2]; % order: ankle, knee, and hip
y_axis_lim=[-45 20;-40 40; -60 30];
x_axis_lim=[-30 20;-15 40; 0 75]; % rajagopal

dataPoint_step   = 1; % just for visualization
marker_list      = ["o","o","o"];
jointLeg_label   = {'ankle','hip','knee'};
description_label= {'knee flexion  0°' 'knee flexion 15°' 'knee flexion 60°';
                   'knee flexion 15°' 'knee flexion 60°' ' ';
                   {'ankle dorsiflexion 20°'; 'hip neutral 0°'} {'ankle plantarflexion 15°'; 'hip neutral 0°'} {'ankle dorsiflexion 20°'; 'hip extension 15°'}};

for trial_opt=1:Misc.nTrials
    trial=trial_list(trial_opt);
    for i=1:trial_sel(trial)
        
        frame_per_trial=50;
        sel_off=(i-1)*frame_per_trial;
        plot_sel=trial_opt+(i-1)*3;
        
        subplot(3,3,plot_sel)
        hold on;
        plot(DatStore(trial).q_exp(1+sel_off:dataPoint_step:frame_per_trial+sel_off),DatStore(trial).T_exp(1+sel_off:dataPoint_step:frame_per_trial+sel_off),marker_list(i),'MarkerEdgeColor',"#000000",'MarkerFaceColor','k','markerSize',5);
        plot(DatStore(trial).q_exp(1+sel_off:dataPoint_step:frame_per_trial+sel_off),genParam(trial).Tsim(1+sel_off:dataPoint_step:frame_per_trial+sel_off),marker_list(i),'MarkerEdgeColor',"#000000",'MarkerFaceColor',"#0072BD",'markerSize',5); hold on;
        plot(DatStore(trial).q_exp(1+sel_off:dataPoint_step:frame_per_trial+sel_off),calParam(trial).Tsim(1+sel_off:dataPoint_step:frame_per_trial+sel_off),marker_list(i),'MarkerEdgeColor',"#000000",'MarkerFaceColor',"#A2142F",'markerSize',5); hold on;
        
        xlabel('angle [deg]');
        ylabel('moment [Nm]')
        set(gca,'FontSize',12);

        if i==1 && (trial_opt==1 || trial_opt==2 || trial_opt==3)
            title(jointLeg_label{trial},'FontSize',20);
        end

        text(x_axis_lim(trial,1)+2,y_axis_lim(trial,1)+20,description_label{trial+(i-1)*3},'fontSize',12)
        ylim(y_axis_lim(trial,:));
        xlim(x_axis_lim(trial,:));
    
    end
end

subplot(3,3,9)
hold on;
plot(0,0,marker_list(i),'MarkerEdgeColor',"#000000",'MarkerFaceColor','k','markerSize',10);
plot(0,0,marker_list(i),'MarkerEdgeColor',"#000000",'MarkerFaceColor',"#A2142F",'markerSize',10); hold on;
plot(0,0,marker_list(i),'MarkerEdgeColor',"#000000",'MarkerFaceColor',"#0072BD",'markerSize',10); hold on;

axis([10,11,10,11]) %move dummy points out of view
legend('experimental','calibrated','generic','location','southwest','fontSize',15); legend boxoff  
axis off %hide axis

function [muscle_coupling_passive]= muscleCouplingPassive(model_conf,dof_name,side)
% Ankle joint
% {'soleus_r'} {'tib_post_r'} {'tib_ant_r'} are separated
% this group is plantar: 'per_brev' and 'per_long', and 'per_tert_r' is
% dorsi, thus per was not coupled
muscle_coupling.ankle_angle.model2392={
    'per_brev' 'per_long';
    'flex_dig' 'flex_hal';
    'ext_dig' 'ext_hal';
    'med_gas' 'lat_gas';
    };
muscle_coupling.ankle_angle.modelrajagopal={
    'perbrev' 'perlong'
    'fdl'     'fhl'
    'edl'     'ehl'
    'gasmed'  'gaslat'
};

% Knee joint
% {'sar_r'}    {'tfl_r'}    {'grac_r'}    {'rect_fem_r'} are separated
% this group is biarticular: 'semimem' 'semiten' 'bifemlh', and bifemsh is
% uniarticular, thus they were not coupled
muscle_coupling.knee_angle.model2392={
    'semimem' 'semiten';
    'semimem' 'bifemlh';
    'vas_med' 'vas_int';
    'vas_med' 'vas_lat';
    };
muscle_coupling.knee_angle.modelrajagopal={
    'semimem' 'semiten';
    'semimem' 'bflh';
    'vasmed'  'vasint';
    'vasmed'  'vaslat';
    };
% Hip joint
% this muscle group is joint: {'add_long_r'} {'add_brev_r'} is separated , and add_mag is too broad
% other muscles have too different muscle paths
muscle_coupling.hip_flexion.model2392={
    'glut_max1' 'glut_max2';
    'glut_max1' 'glut_max3';
    'glut_med1' 'glut_med2';
    'glut_med1' 'glut_med3';
    'glut_min1' 'glut_min2';
    'glut_min1' 'glut_min3';
    'add_mag1'  'add_mag2';
    'add_mag1'  'add_mag3';
    };
muscle_coupling.hip_flexion.modelrajagopal={
    'glmax1' 'glmax2';
    'glmax1' 'glmax3';
    'glmed1' 'glmed2';
    'glmed1' 'glmed3';
    'glmin1' 'glmin2';
    'glmin1' 'glmin3';
    'addmagDist'  'addmagIsch';
    'addmagDist'  'addmagMid';
    'addmagDist'  'addmagProx'; % one more than in 2392
    };

muscle_coupling_passive=append(muscle_coupling.(dof_name).(['model' model_conf{:}]),['_' side]); 
end