function [Misc] = getMuscleProperties(model_path,Misc)
% --------------------------------------------------------------------------
%getMuscleProperties
%     Reads properties of the muscle from the .osim file and handles user
%     defined tendon stiffness
% 
% INPUT:
%     model_path
%     Path to the .osim model
% 
%     Misc
%     Miscellaneous info used through the code
% 
% OUTPUT:
%     Misc
%     Miscellaneous info used through the code
% 
% Original author: 
% Original date: 
%
% Last edit by: Dhruv Gupta
% Last edit date: May 3, 2022
% --------------------------------------------------------------------------

% read the model
import org.opensim.modeling.*;
model = Model(model_path);

% read the muscle properties
muscles = model.getMuscles();
nNames = muscles.getSize();
params = zeros(8, nNames);

PassiveParam_kpe=Misc.PassiveParam_kpe;
PassiveParam_ksf=Misc.PassiveParam_ksf;
PassiveParam_ke0=Misc.PassiveParam_ke0;

for i = 1:nNames
   muscle = muscles.get(i-1);
   Misc.allMuscleList{i} = char(muscle.getName());
   params(3,i) = muscle.getTendonSlackLength();		
   params(2,i) = muscle.getOptimalFiberLength(); 	
   params(1,i) = muscle.getMaxIsometricForce();  	
   params(4,i) = muscle.getPennationAngleAtOptimalFiberLength(); 
   params(5,i) = muscle.getMaxContractionVelocity();
   params(6,i) = PassiveParam_kpe(i); % kpe default value 4
   params(7,i) = PassiveParam_ksf(i); % ksf default value 1
   params(8,i) = PassiveParam_ke0(i); % ke0 default value 0.6
end

Misc.nAllMuscList = length(Misc.allMuscleList);

% create additional variables with the same information
Misc.FMo=params(1,:);
Misc.lMo=params(2,:);
Misc.lTs=params(3,:);
Misc.alphao=params(4,:);
Misc.kpe=params(6,:);
Misc.ksf=params(7,:);
Misc.ke0=params(8,:);
Misc.params = params;


if ~isfield(Misc,'kT') || isempty(Misc.kT)
    Misc.kT =ones(1,length(Misc.allMuscleList)).*35;
end

% set the default value of the tendon stiffness
if isfield(Misc,'Set_kT_ByName') && ~isempty(Misc.Set_kT_ByName)
    Misc = set_kT_ByName(Misc);
end

% Shift tendon force-length curve as a function of the tendon stiffness
Misc.shift = getShift(Misc.kT);

end