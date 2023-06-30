function[MC_parameter,E_rate_SCALED,w_rate,h_rate_SCALED]= MC_LW05_R(muscle_parameter,muscle_DynCon,option,basalOn)
% read inputs
f_FT= muscle_parameter.FT;
OFL = muscle_parameter.OFL;
MCV = muscle_parameter.MCV;
MIF = muscle_parameter.MIF;                  
PCSA= muscle_parameter.PCSA;
MASS= muscle_parameter.mass;

excitation= muscle_DynCon.muscle_excitation;
activation= muscle_DynCon.muscle_activation;
MTUforce  = muscle_DynCon.muscle_MTUforce;    
V_CE      = muscle_DynCon.V_CE; %+shortening -lengthening
W_CE      = muscle_DynCon.W_CE; % W_CE= F_CE*V_CE: +W_CE spent -W_CE gain

lMtilde  = muscle_DynCon.lMtilde;
vMtilde  = muscle_DynCon.vMtilde;
fl_act_multiplier=muscle_DynCon.fl_act_multiplier;
f_v_multiplier   =muscle_DynCon.f_v_multiplier;
fl_pas_multiplier=muscle_DynCon.fl_pas_multiplier;

time=muscle_DynCon.time;
data_length= length(time);
%% STABLE HEAT RATE
% Constants
G      = 4; % curvature of the force-velocity curve

% factors
norVal  = MIF*OFL; % N*m
P= MIF.*activation.*fl_act_multiplier.*f_v_multiplier;

h_M_rate       = activation.*MCV./(G^2); % OFL/s
h_M_rate_SCALED= h_M_rate*norVal; % N*m/s

c=0; 
t_stim=double.empty;
for i=1:data_length
     if excitation(i)>0.1
        c=c+1;
        t_stim(i,1)=c; 
     else
        t_stim(i,1)=0;
        c=0;
     end
end
time_dc_unit=(time(end)-time(1))/(length(time)-1);
t_stim=t_stim*time_dc_unit;

%% LABILE HEAT RATE
h_L_rate= h_M_rate.*(0.8.*exp(-0.72.*t_stim(:))+0.175.*exp(-0.022.*t_stim(:))); % OFL/s
h_L_rate_SCALED = h_L_rate.*norVal; % N*m/s
% h_L_rate_SCALED= zeros(data_length,1);

%% SHORTENING/LENGTHENING HEAT RATE
h_SL_rate = zeros(data_length,1);
for i=1:data_length
    if V_CE(i) >=0 % shortening
        h_SL_rate(i)= activation(i).*(V_CE(i)/OFL)./G; % OFL/s
    end
end
h_SL_rate_SCALED= h_SL_rate*norVal; % N*m/s

%% THERMOELASTIC HEAT RATE
muscle_force_smooth = smooth(time,P,0.05,'loess');
muscle_force_rate = gradient(muscle_force_smooth)./gradient(time);
muscle_force_rate(isnan(muscle_force_rate))=0;
h_T_rate= -0.014*muscle_force_rate; % N/s absortion from the muscle
h_T_rate_SCALED= h_T_rate*OFL;

%% WORK RATE
w_rate = W_CE;

%% ENERGY RATE
h_rate_SCALED=zeros(data_length,1);
E_rate_SCALED=zeros(data_length,1);
for i=1:data_length
    if V_CE(i)>=0 %shortening
         h_rate_SCALED(i,1)= h_M_rate_SCALED(i) + h_L_rate_SCALED(i) + h_SL_rate_SCALED(i) + h_T_rate_SCALED(i);
         E_rate_SCALED(i,1)= w_rate(i) + h_rate_SCALED(i,1);
    
    elseif V_CE(i)<=0 %lengthening
         h_rate_SCALED(i,1)= h_M_rate_SCALED(i).*(0.3+0.7.*exp(-8*((P(i)./activation(i))-1) ))+w_rate(i);
         E_rate_SCALED(i,1)= w_rate(i) + h_rate_SCALED(i);
        
    % h_rate_SCALED(i)= h_M_rate_SCALED(i).*(0.3+0.7.*exp(-8*((P(i)./ACT(i))-1) ))-w_rate(i);
    % I assigned the scaling factor to the mantainance term (replacing the previous equation). 
    % All the others are zero. I do this only because I have to assign this to a
    % particular heat. It is not written in the original formula but it
    % makes sense
    
        h_M_rate_SCALED(i,1) = h_M_rate_SCALED(i).*(0.3+0.7.*exp(-8*((P(i)./activation(i))-1) ))+w_rate(i);
        h_SL_rate_SCALED(i,1)=0;
        h_T_rate_SCALED(i,1)=0;
        h_L_rate_SCALED(i,1)=0;
    end
end

E_value= cumtrapz(time,E_rate_SCALED);
h_rate_AM=h_M_rate_SCALED + h_T_rate_SCALED + h_L_rate_SCALED;
h_rate_SL=h_SL_rate_SCALED;

MC_parameter(1,:)= E_rate_SCALED;
MC_parameter(2,:)= w_rate;
MC_parameter(3,:)= h_rate_SCALED;
MC_parameter(4,:)= h_rate_AM;
MC_parameter(5,:)= h_rate_SL;

MC_parameter(6,:)= h_M_rate_SCALED;
MC_parameter(7,:)= h_SL_rate_SCALED;
MC_parameter(8,:)= h_T_rate_SCALED;
MC_parameter(9,:)= E_value;
end