function[MC_parameter,E_rate,w_rate,h_rate]= MC_UC16_R(muscle_parameter,muscle_DynCon,option,basalOn)
% read inputs
f_FT= muscle_parameter.FT;      % from 0 to 1
OFL = muscle_parameter.OFL;
MCV = muscle_parameter.MCV;
MIF = muscle_parameter.MIF;                 
PCSA= muscle_parameter.PCSA;
MASS= muscle_parameter.mass;

excitation= muscle_DynCon.muscle_excitation;
activation= muscle_DynCon.muscle_activation;
MTUforce  = muscle_DynCon.muscle_MTUforce;     
V_CE      = muscle_DynCon.V_CE; %-shortening +lengthening
W_CE      = muscle_DynCon.W_CE; % W_CE=-F_CE*V_CE: +W_CE spent -W_CE gain

lMtilde  = muscle_DynCon.lMtilde;
vMtilde  = muscle_DynCon.vMtilde;
fl_act_multiplier=muscle_DynCon.fl_act_multiplier(:);
f_v_multiplier   =muscle_DynCon.f_v_multiplier(:);
fl_pas_multiplier=muscle_DynCon.fl_pas_multiplier(:);

time=muscle_DynCon.time;
data_length= length(time);
%% Physiological related parameters
% rho= 1059.7;    sigma=0.25*10^6;
% m= PCSA*rho*muscle_OFL;
S= 1.5; %Aerobic

VCEtilde= V_CE./OFL; % [OFL/s]
L_CE= lMtilde*OFL;   % [m]
f_FT  = muscle_parameter.FT.*ones(data_length,1);
f_ST  = (1-f_FT);
%% SCALING FACTOR
A=zeros(data_length,1);
for i=1:data_length
    if excitation(i) > activation(i)
    A(i,1)= excitation(i);
    elseif excitation(i) <= activation(i)
    A(i,1)= (excitation(i)+activation(i))/2;      
    end
end
A_AM= A.^(0.6);
A_S = zeros(data_length,1);
for i=1:data_length
    if VCEtilde(i) <= 0 %shortening -
    A_S(i,1)= A(i).^2;
    elseif VCEtilde(i) > 0 %lengthening +
    A_S(i,1)= A(i);
    end
end 
%% VELOCITY PARAMETERS
VCE_MAX= MCV;

VCE_MAX_ST= VCE_MAX/2.5; %4
VCE_MAX_FT= VCE_MAX; %10

alpha_S_ST= 100/VCE_MAX_ST;
alpha_S_FT= 153/VCE_MAX_FT;
alpha_L= 4*alpha_S_ST; % umberger 2003
%% MUSCLE RECRUITMENT
% Equations
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

f_active_FT=1-f_active_ST;%
%% ACTIVATION/MAINTENANCE HEAT RATE
h_AM_rate_scaled=zeros(data_length,1);
h_AM_rate= 128.*f_active_FT + 25; %isometric

for i=1:data_length
    if L_CE(i)<=OFL
        h_AM_rate_scaled(i,1)= h_AM_rate(i).*A_AM(i).*S*MASS;
    elseif L_CE(i)>OFL
        h_AM_rate_scaled(i,1)= (0.4*h_AM_rate(i)+0.6*h_AM_rate(i)*fl_act_multiplier(i)).*A_AM(i).*S*MASS;
    end
end
%% Shortening/lengthening heat rate
h_SL_rate_scaled=zeros(data_length,1);
h_SL_rate=zeros(data_length,1);

for i=1:data_length
   if VCEtilde(i)<=0 % shortening -
      h_SL_rate(i)= -(alpha_S_ST.*VCEtilde(i).*f_active_ST(i) + alpha_S_FT.*VCEtilde(i).*f_active_FT(i)).*A_S(i);
   elseif VCEtilde(i)>0 % lengthening +
      h_SL_rate(i)= alpha_L.*VCEtilde(i).*A(i);
   end
end

for i=1:data_length
    if L_CE(i)<=OFL
        h_SL_rate_scaled(i,1)= h_SL_rate(i)                      *S*MASS;
    elseif L_CE(i)>OFL
        h_SL_rate_scaled(i,1)= h_SL_rate(i).*fl_act_multiplier(i)*S*MASS;
    end
end
%% WORK RATE
w_rate = W_CE;
%% ENERGY RATE
h_SL_rate_scaled_updated=zeros(data_length,1);

E_rate_previous= h_AM_rate_scaled+h_SL_rate_scaled+w_rate;
for i=1:data_length % net negative energy is constraint/limited
    if E_rate_previous(i)<=0
    h_SL_rate_scaled_updated(i,1)= -h_AM_rate_scaled(i)-w_rate(i);
    else
    h_SL_rate_scaled_updated(i,1)= h_SL_rate_scaled(i);   
    end
end

h_rate_scaled= h_AM_rate_scaled+h_SL_rate_scaled_updated;
E_rate       = w_rate+h_rate_scaled;

if basalOn==1
    h_rate_scaled=h_rate_scaled+1.2*MASS; %1.2 W/kg
    E_rate=h_rate_scaled + w_rate;
end

h_AM_rate=h_AM_rate_scaled;
h_SL_rate=h_SL_rate_scaled_updated;
h_rate   =h_rate_scaled;

if data_length>1
    E_value= cumtrapz(time,E_rate);
else
    E_value=nan;
end

h_rate_AM=h_AM_rate_scaled;
h_rate_SL=h_SL_rate_scaled_updated;

MC_parameter(1,:)= E_rate;
MC_parameter(2,:)= w_rate;
MC_parameter(3,:)= h_rate;
MC_parameter(4,:)= h_rate_AM;
MC_parameter(5,:)= h_rate_SL;

MC_parameter(6,:)= h_AM_rate;
MC_parameter(7,:)= h_SL_rate;
MC_parameter(8,:)= E_value;
end