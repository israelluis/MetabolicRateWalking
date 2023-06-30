function[MC_parameter,E_rate_SCALED,w_rate,h_rate_SCALED]= MC_LW07_R(muscle_parameter,muscle_DynCon,option,basalOn)
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
V_CE=muscle_DynCon.V_CE; %-shortening +lengthening
W_CE=muscle_DynCon.W_CE; % W_CE=-F_CE*V_CE where W_CE +GENERATION(spending) AND -ABSORTION(gaining)

lMtilde  = muscle_DynCon.lMtilde;
vMtilde  = muscle_DynCon.vMtilde;
fl_act_multiplier=muscle_DynCon.fl_act_multiplier;
f_v_multiplier   =muscle_DynCon.f_v_multiplier;
fl_pas_multiplier=muscle_DynCon.fl_pas_multiplier;

time=muscle_DynCon.time;
data_length= length(time);
%% ACTIVATION/MAINTENANCE HEAT RATE
% MAINTENANCE
% Constants
G      = 4; % curvature of the force-velocity curve
gamma  = 1.5; % According to LW07 SI

% factors
X      = activation.*fl_act_multiplier; % scaling factor
norVal = MIF*OFL; % unnormalizing factor Po*Lo [N*m]
P      = f_v_multiplier;

h_M_rate=zeros(data_length,1);

for i=1:data_length
%     if V_CE(i) >0 % lengthening
     if f_v_multiplier(i)>=1 % lengthening
        h_M_rate(i,1)=0.3*gamma*MCV./(G.^2)+0.7*gamma*MCV./(G.^2)*exp(-7*MCV*(P(i)-1)); % OFL/s (positive)
    else
        h_M_rate(i,1)=gamma*MCV./(G.^2); % OFL/s
    end
end
h_M_rate_scaled_heat= (0.3.*activation.*h_M_rate+0.7*X.*h_M_rate).*norVal; % (1/s)* N*m =N*m/s

%% SHORTENING/LENGTHENING HEAT RATE
h_SL_rate =zeros(data_length,1);
for i=1:data_length
    if V_CE(i) >0 % lengthening +
        h_SL_rate(i,1)= 0.5*P(i).*V_CE(i)/OFL; % 1/s (positive) (it becomes 50%W_CE dissipated as heat)
    elseif V_CE(i) <=0 % shortening -
        h_SL_rate(i,1)=-(V_CE(i)./OFL)./G; % 1/s
    end
end
h_SL_rate_scaled_heat= h_SL_rate.*X.*norVal; % 1/s * N*m = N*m/s

%% WORK RATE
w_rate = W_CE;

%% ENERGY RATE
h_rate_SCALED= h_M_rate_scaled_heat+h_SL_rate_scaled_heat;
E_rate_SCALED= w_rate+h_rate_SCALED;

if basalOn==1
    h_rate_SCALED=h_rate_SCALED+1.2*MASS; %1.2 W/kg
    E_rate_SCALED=h_rate_SCALED + w_rate;
end
h_rate_AM=h_M_rate_scaled_heat;
h_rate_SL=h_SL_rate_scaled_heat;

if data_length>1
    E_value  =cumtrapz(time,E_rate_SCALED);
else
    E_value=nan;
end


MC_parameter(1,:)= E_rate_SCALED;
MC_parameter(2,:)= w_rate;
MC_parameter(3,:)= h_rate_SCALED;
MC_parameter(4,:)= h_rate_AM;
MC_parameter(5,:)= h_rate_SL;

MC_parameter(6,:)= h_M_rate_scaled_heat;
MC_parameter(7,:)= h_SL_rate_scaled_heat;
MC_parameter(8,:)= E_value;
end