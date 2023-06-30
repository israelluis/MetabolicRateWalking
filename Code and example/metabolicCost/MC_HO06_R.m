function[MC_parameter,E_rate,w_rate,h_rate]= MC_HO06_R(muscle_parameter,muscle_DynCon,option,basalOn)
% read inputs
f_FT= muscle_parameter.FT;
OFL = muscle_parameter.OFL;
MCV = muscle_parameter.MCV;
MIF = muscle_parameter.MIF;                      % NEEDED
PCSA= muscle_parameter.PCSA;
MASS= muscle_parameter.mass;

excitation= muscle_DynCon.muscle_excitation;
activation= muscle_DynCon.muscle_activation;
MTUforce  = muscle_DynCon.muscle_MTUforce;       % NOT NEEDED
V_CE      = muscle_DynCon.V_CE; %-shortening +lengthening
W_CE      = muscle_DynCon.W_CE; % W_CE=-F_CE*V_CE: +W_CE spent -W_CE gain

lMtilde  = muscle_DynCon.lMtilde;
vMtilde  = muscle_DynCon.vMtilde;
fl_act_multiplier=muscle_DynCon.fl_act_multiplier;
f_v_multiplier   =muscle_DynCon.f_v_multiplier;
fl_pas_multiplier=muscle_DynCon.fl_pas_multiplier;

time=muscle_DynCon.time;
data_length= length(time);
%% ANNOTATION - You need Hatse and buys 1977 to implement this
% we know that fi_hatze=gDotTilde/(gDotTilde+hDotTilde) and also
% fi_hatze = for FT -> 0.35 w/kg and for ST -> 0.45 w/kg
% gDotTilde+hDotTilde = for FT -> 150 w/kg and for ST -> 24.4 w/kg 
% Thus gDotTilde = for FT -> 52.50 w/kg and for ST -> 10.98 w/kg 
f_FT  = muscle_parameter.FT.*ones(data_length,1);
f_ST  = (1-f_FT);
f_active_ST = f_ST;
f_active_FT = f_FT;
%% MUSCLE RECRUITMENT FROM Uchida et al. 2016
if option==2
u_f   = 1-cos(excitation*pi/2);
u_s   = sin(excitation*pi/2);
f_active_ST = zeros(data_length,1); % 0-1
for i=1:data_length
    if excitation(i)==0
    f_active_ST(i,1)=1;
    else
    f_active_ST(i,1)= f_ST(i)*u_s(i)./(f_ST(i).*u_s(i)+f_FT(i).*u_f(i));   
    end
end
f_active_FT=1-f_active_ST;
end
%% ACTIVATION/MAINTENANCE HEAT RATE
% ACTIVATION HEAT RATE
% Constants
h_A_ST=10.98; % w/kg
h_A_FT=52.50; % w/kg
k3_FT= 12;  k4_FT= 14; % from hatze (Table1)
k3_ST= 6;   k4_ST= 8;  % from hatze (Table1)

h_A= h_A_FT.*f_active_FT+h_A_ST.*f_active_ST;
k3= k3_FT.*f_active_FT+k3_ST.*f_active_ST;
k4= k4_FT.*f_active_FT+k4_ST.*f_active_ST;

vMAXFreq= k3+k4.*activation; % hazte equation 20
vRelativeStimulationFreq= activation.^2; % activation*(stimulation_rate/max_stimulation_rate)=activation*activation=activation.^2

h_A_rate= MASS*h_A.*vRelativeStimulationFreq.*( (1-exp(-0.25-18.2./(vRelativeStimulationFreq.*vMAXFreq)))./(1-exp(-0.25-18.2./vMAXFreq)) ); 
%% MAINTANCE HEAT RATE
% Constants
h_M_ST=13.42; % w/kg % from h_A_f and h_A+h_M
h_M_FT=97.50; % w/kg % from h_A_f and h_A+h_M

h_M= h_M_FT.*f_active_FT+h_M_ST.*f_active_ST;
h_A_f= h_A./(h_A+h_M);

h_M_rate= MASS.*(h_A+h_M).*activation.*(fl_act_multiplier-h_A_f);
%% SHORTENING/LENGTHENING HEAT RATE
aTilde_coefficient_FT=0.28; % w/kg
aTilde_coefficient_ST=0.16; % w/kg

aTilde= (aTilde_coefficient_FT.*f_active_FT+aTilde_coefficient_ST.*f_active_ST).*MIF;

h_SL_rate=zeros(data_length,1);
for i=1:data_length
   if V_CE(i)<=0 %shortening -
      h_SL_rate(i,1)= -aTilde(i).*activation(i).*fl_act_multiplier(i).*V_CE(i); % shortening depend approximately linearly on absolute shortening velocity
   elseif V_CE(i)>0 %lengthening +
      h_SL_rate(i,1)= 0; % there is no definition of lenghtening in HO06
   end
end
%% WORK RATE
w_rate = W_CE;
%% ENERGY RATE
h_rate= h_A_rate + h_M_rate + h_SL_rate;
E_rate= w_rate + h_rate;
    
if option==0
    % keep it as it was :D
elseif option>=1 
    for i=1:data_length
    if E_rate(i)<=0
       h_SL_rate(i,1)= -h_A_rate(i) - h_M_rate(i) - w_rate(i);
    end
    end
    h_rate= h_A_rate + h_M_rate + h_SL_rate;
    E_rate= h_rate + w_rate;
end

if basalOn==1
    h_rate=h_rate+1.2*MASS; %1.2 W/kg
    E_rate=h_rate + w_rate;
end

h_rate_AM=h_A_rate + h_M_rate;
h_rate_SL=h_SL_rate;

if data_length>1
    E_value  =cumtrapz(time,E_rate);
else
    E_value=nan;
end

MC_parameter(1,:)= E_rate;
MC_parameter(2,:)= w_rate;
MC_parameter(3,:)= h_rate;
MC_parameter(4,:)= h_rate_AM;
MC_parameter(5,:)= h_rate_SL;

MC_parameter(6,:)= h_A_rate;
MC_parameter(7,:)= h_M_rate;
MC_parameter(8,:)= h_SL_rate;
MC_parameter(9,:)= E_value;
end