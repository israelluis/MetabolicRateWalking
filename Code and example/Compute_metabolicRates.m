%% METABOLIC RATES
%
% This script serves as an example to run the computation of metabolic rates. 
%% DESCRIPTION
% 
% We first setup the information about the subject, and the results
% obtained by any of the simulation workflows. Then we organize the 
% results and compute the metabolic rates based on metabolic energy models.
% Finally we plot the metabolic rate, work rate, and heat rate.
%
%% SETUP SUBJECT INFO AND READ SIMULATION WORKFLOW
clc;
MainPath = pwd;

% Subject data
subjectMass  = 72.6;
subjectHeight=174.3;
    
% Compute leg volumen based on handsfield 2014
V_total_leg    = 47*subjectMass*(subjectHeight/100)+1285; % Input: Kg & meter - Output: cm^3
rho            = 1059.7;                                   % kg/m^3 [Umberger et al. 2003]

% Load result
simulation_workflow='effort_gen';   % choose simulation workflow: effort_gen effort_cal emg_ten
MRS    =load(fullfile(MainPath,'Results',simulation_workflow,'simulation_Results.mat')); % add the simulation workflow path
Misc   =MRS.Misc;
Results=MRS.Results;

% Get simulation workflow label from the MRS
simulation_type= char(fieldnames(Results.Time));
NTrials        = length(Misc.MAtrialName);
NMuscles       = length(Results.MuscleNames);
%% COMPUTE METABOLIC RATES

mode_basalOn=1; % Boolean: Select if you want to add a basal heat
mode_negWork=1; % Boolean: Select if you want negative heat to be dissipated, thus net energy rate is >=0
        
fr_ext=5; % take out the first and last 5 frames, MRS does not guarantee optimality at the beginning and end of gait cycle
MetCost=cell(NTrials,NMuscles);

for trial_sel=1:NTrials
    for muscle_sel=1:NMuscles
        
        % organize muscle-tendon parameters
        muscle_name    = string(Results.MuscleNames);

        muscle_parameter.MIF = Misc.params(1,muscle_sel)';
        muscle_parameter.OFL = Misc.params(2,muscle_sel)';
        muscle_parameter.MCV = Misc.params(5,muscle_sel)';

        muscle_parameter.VOL = getVolumeFraction(muscle_name(muscle_sel))/100*(V_total_leg./(100)^3);    % m^3 
        muscle_parameter.PCSA= (muscle_parameter.VOL)./(muscle_parameter.OFL);
        muscle_parameter.FT  = (1-getSlowTwitchRatios_Upd(muscle_name(muscle_sel))); % [fraction]
        muscle_parameter.mass= muscle_parameter.VOL.*rho;

        % organize muscle-tendon states
        iSel = 1+fr_ext:length(Results.MActivation(trial_sel).(simulation_type)(1,:))-1-fr_ext; % to take out the first and final frames

        muscle_DynCon.time             = Results.Time(trial_sel).(simulation_type)(iSel,:); 
        muscle_DynCon.muscle_excitation= Results.MExcitation(trial_sel).(simulation_type)(muscle_sel,iSel)'; 
        muscle_DynCon.muscle_activation= Results.MActivation(trial_sel).(simulation_type)(muscle_sel,iSel)';

        muscle_DynCon.lMtilde          = Results.lMtildeopt(trial_sel).(simulation_type)(muscle_sel,iSel)';
        muscle_DynCon.vMtilde          = (Results.vMtilde(trial_sel).(simulation_type)(muscle_sel,iSel)./muscle_parameter.MCV)';

        muscle_DynCon.fl_act_multiplier = Results.FMltilde(trial_sel).(simulation_type)(muscle_sel,iSel)'; %FMltilde []
        muscle_DynCon.f_v_multiplier    = Results.FMvtilde(trial_sel).(simulation_type)(muscle_sel,iSel)'; %FMvtilde []

        muscle_DynCon.muscle_MTUforce  = Results.TForce(trial_sel).(simulation_type)(muscle_sel,iSel)';
        muscle_DynCon.fl_pas_multiplier= (Results.Fpe(trial_sel).(simulation_type)(muscle_sel,iSel)./muscle_parameter.MIF)'; %Fpe [N]

        muscle_DynCon.F_CE = (muscle_DynCon.muscle_activation).*(muscle_DynCon.fl_act_multiplier).*(muscle_DynCon.f_v_multiplier).*muscle_parameter.MIF;
        muscle_DynCon.V_CE = (muscle_DynCon.vMtilde.*muscle_parameter.OFL.*muscle_parameter.MCV);
        muscle_DynCon.W_CE =-(muscle_DynCon.F_CE.*muscle_DynCon.V_CE);

        % metabolic cost models
        [MC_parameter_UM03_0,~,~,~]= MC_UM03_R(muscle_parameter,muscle_DynCon,mode_negWork,mode_basalOn);
        MetCost(trial_sel,muscle_sel,1)= {MC_parameter_UM03_0};
        
        [MC_parameter_BH04_0,~,~,~]  = MC_BH04_R(muscle_parameter,muscle_DynCon,mode_negWork,mode_basalOn);
        MetCost(trial_sel,muscle_sel,2)= {MC_parameter_BH04_0};
        
        [MC_parameter_HO06_1,~,~,~]= MC_HO06_R(muscle_parameter,muscle_DynCon,mode_negWork,mode_basalOn);
        MetCost(trial_sel,muscle_sel,3)= {MC_parameter_HO06_1};
        
        [MC_parameter_LW07_0,~,~,~]= MC_LW07_R(muscle_parameter,muscle_DynCon,mode_negWork,mode_basalOn);
        MetCost(trial_sel,muscle_sel,4)= {MC_parameter_LW07_0};
        
        [MC_parameter_UM10_0,~,~,~]= MC_UM10_R(muscle_parameter,muscle_DynCon,mode_negWork,mode_basalOn);
        MetCost(trial_sel,muscle_sel,5)= {MC_parameter_UM10_0};
        
        [MC_parameter_UC16_0,~,~,~]= MC_UC16_R(muscle_parameter,muscle_DynCon,mode_negWork,mode_basalOn);
        MetCost(trial_sel,muscle_sel,6)= {MC_parameter_UC16_0};
    end
end
%% PLOT MUSCLE METABOLIC RATES
fig=figure(1); clf;
set(gcf,'color','w','Visible','on','Position',[50 50 1200 700]);
label_metCost={'Umberger et al. 2003' 'Bhargava et al. 2004' 'Houdijk et al. 2006' ...
               'Lichtwark & Wilson 2007' 'Umberger 2010' 'Uchida et al. 2016'};
trial_sel  =2; % walking speed: 1 slow, 2 normal, 3 fast
metCost_sel=2; % choose the metabolic energy model

for muscle_sel=1:NMuscles
    subplot(5,8,muscle_sel)
    hold on
    muscle_energetics=MetCost{trial_sel,muscle_sel,metCost_sel};
    data_length=length(muscle_energetics(1,:));
    gait_cycle=linspace(0,100,data_length); % just for visualization
    plot(gait_cycle,muscle_energetics(1,:),'k','displayName','Met'); % metabolic rate
    plot(gait_cycle,muscle_energetics(2,:),'b','displayName','Work'); % work rate
    plot(gait_cycle,muscle_energetics(3,:),'r','displayName','Heat'); % heat rate
    title(Results.MuscleNames{muscle_sel});
    xlabel('gait cycle');
    ylabel('rate [J]');
end
sgtitle(['Simulation Workflow: ' simulation_workflow ' - Metabolic Cost: ' label_metCost{metCost_sel}],'interpreter','none');