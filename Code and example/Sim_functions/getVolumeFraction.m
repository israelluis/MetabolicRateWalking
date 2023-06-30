% This function returns the volume fraction in the muscles
% The data come from Handsfield (2014). Provide % of the fraction
%
% Author: Israel Luis
% Date: 27/01/2021
function mus_vol = getVolumeFraction(muscleNames) 

    sigma_data.glut_med1_r = 4.54/3;   sigma_data.glmed1_r = 4.54/3;
    sigma_data.glut_med2_r = 4.54/3;   sigma_data.glmed2_r = 4.54/3;
    sigma_data.glut_med3_r = 4.54/3;   sigma_data.glmed3_r = 4.54/3;
    sigma_data.glut_min1_r = 1.47/3;    sigma_data.glmin1_r = 1.47/3;
    sigma_data.glut_min2_r = 1.47/3;    sigma_data.glmin2_r = 1.47/3;
    sigma_data.glut_min3_r = 1.47/3;    sigma_data.glmin3_r = 1.47/3;
    sigma_data.semimem_r = 3.46;     
    sigma_data.semiten_r = 2.60;    
    sigma_data.bifemlh_r = 2.92;     sigma_data.bflh_r = 2.92;        sigma_data.bflh140_r = 2.92;
    sigma_data.bifemsh_r = 1.40;     sigma_data.bfsh_r = 1.40;        sigma_data.bfsh140_r = 1.40; 
    sigma_data.sar_r = 2.29;         sigma_data.sart_r = 2.29;
    sigma_data.add_mag1_r = 7.86/3;    sigma_data.addmagDist_r = 7.86/3;
    sigma_data.add_mag2_r = 7.86/3;    sigma_data.addmagIsch_r = 7.86/3;
    sigma_data.add_mag3_r = 7.86/3;    sigma_data.addmagMid_r = 7.86/3;     sigma_data.addmagProx_r = 7.86/3;
    sigma_data.tfl_r = 0.89;
    sigma_data.pect_r = 0.92;
    sigma_data.grac_r = 1.46;
    sigma_data.glut_max1_r = 11.93/3;   sigma_data.glmax1_r = 11.93/3;
    sigma_data.glut_max2_r = 11.93/3;   sigma_data.glmax2_r = 11.93/3;
    sigma_data.glut_max3_r = 11.93/3;   sigma_data.glmax3_r = 11.93/3;
    sigma_data.iliacus_r = 2.48;
    sigma_data.psoas_r = 3.80;       
    sigma_data.quad_fem_r = 0.45;
    sigma_data.gem_r = 0.23; % gemelli not described
    sigma_data.peri_r = 0.61;         sigma_data.piri_r = 0.61;
    sigma_data.rect_fem_r = 3.79;     sigma_data.recfem_r = 3.79;
    sigma_data.vas_med_r = 6.06;      sigma_data.vasmed_r = 6.06;
    sigma_data.vas_int_r = 3.84;      sigma_data.vasint_r = 3.84;
    sigma_data.vas_lat_r = 11.66;     sigma_data.vaslat_r = 11.66;      sigma_data.vaslat140_r = 11.66;
    sigma_data.med_gas_r = 3.62;      sigma_data.gasmed_r = 3.62;
    sigma_data.lat_gas_r = 2.11;      sigma_data.gaslat_r = 2.11;      sigma_data.gaslat140_r = 2.11;
    sigma_data.soleus_r = 6.21;
    sigma_data.tib_post_r = 1.49;     sigma_data.tibpost_r = 1.49;
    sigma_data.flex_dig_r = 0.43;     sigma_data.fdl_r = 0.43;
    sigma_data.flex_hal_r = 1.09;     sigma_data.fhl_r = 1.09;
    sigma_data.tib_ant_r = 1.91;      sigma_data.tibant_r = 1.91; 
    sigma_data.per_brev_r = 1.83/2;   sigma_data.perbrev_r = 1.83/2; % peroneals (brev/long) are combined
    sigma_data.per_long_r = 1.83/2;   sigma_data.perlong_r = 1.83/2; % peroneals (brev/long) are combined
    sigma_data.per_tert_r = 1.83/3;   % peroneals (tertius) assumed to be 1/3
    sigma_data.ext_dig_r = 0.43;     sigma_data.edl_r = 0.43;
    sigma_data.ext_hal_r = 1.09;     sigma_data.ehl_r = 1.09;
    sigma_data.ercspn_r = 99.0; %Back muscle not described
    sigma_data.intobl_r = 99.0; %Back muscle not described
    sigma_data.extobl_r = 99.0; %Back muscle not described    
    sigma_data.add_long_r = 2.26;    sigma_data.addlong_r = 2.26;
    sigma_data.add_brev_r = 1.47;    sigma_data.addbrev_r = 1.47;
    
    sigma_data.glut_med1_l = 4.54/3;   sigma_data.glmed1_l = 4.54/3;
    sigma_data.glut_med2_l = 4.54/3;   sigma_data.glmed2_l = 4.54/3;
    sigma_data.glut_med3_l = 4.54/3;   sigma_data.glmed3_l = 4.54/3;
    sigma_data.glut_min1_l = 1.47/3;    sigma_data.glmin1_l = 1.47/3;
    sigma_data.glut_min2_l = 1.47/3;    sigma_data.glmin2_l = 1.47/3;
    sigma_data.glut_min3_l = 1.47/3;    sigma_data.glmin3_l = 1.47/3;
    sigma_data.semimem_l = 3.46;     
    sigma_data.semiten_l = 2.60;    
    sigma_data.bifemlh_l = 2.92;     sigma_data.bflh_l = 2.92;        sigma_data.bflh140_l = 2.92;
    sigma_data.bifemsh_l = 1.40;     sigma_data.bfsh_l = 1.40;        sigma_data.bfsh140_l = 1.40; 
    sigma_data.sar_l = 2.29;         sigma_data.sart_l = 2.29;
    sigma_data.add_mag1_l = 7.86/3;    sigma_data.addmagDist_l = 7.86/3;
    sigma_data.add_mag2_l = 7.86/3;    sigma_data.addmagIsch_l = 7.86/3;
    sigma_data.add_mag3_l = 7.86/3;    sigma_data.addmagMid_l = 7.86/3;     sigma_data.addmagProx_l = 7.86/3;
    sigma_data.tfl_l = 0.89;
    sigma_data.pect_l = 0.92;
    sigma_data.grac_l = 1.46;
    sigma_data.glut_max1_l = 11.93/3;   sigma_data.glmax1_l = 11.93/3;
    sigma_data.glut_max2_l = 11.93/3;   sigma_data.glmax2_l = 11.93/3;
    sigma_data.glut_max3_l = 11.93/3;   sigma_data.glmax3_l = 11.93/3;
    sigma_data.iliacus_l = 2.48;
    sigma_data.psoas_l = 3.80;       
    sigma_data.quad_fem_l = 0.45;
    sigma_data.gem_l = 0.23; % gemelli not described
    sigma_data.peri_l = 0.61;         sigma_data.piri_l = 0.61;
    sigma_data.rect_fem_l = 3.79;     sigma_data.recfem_l = 3.79;
    sigma_data.vas_med_l = 6.06;      sigma_data.vasmed_l = 6.06;
    sigma_data.vas_int_l = 3.84;      sigma_data.vasint_l = 3.84;
    sigma_data.vas_lat_l = 11.66;     sigma_data.vaslat_l = 11.66;     sigma_data.vaslat140_l = 11.66;
    sigma_data.med_gas_l = 3.62;      sigma_data.gasmed_l = 3.62;
    sigma_data.lat_gas_l = 2.11;      sigma_data.gaslat_l = 2.11;      sigma_data.gaslat140_l = 2.11;
    sigma_data.soleus_l = 6.21;
    sigma_data.tib_post_l = 1.49;     sigma_data.tibpost_l = 1.49;
    sigma_data.flex_dig_l = 0.43;     sigma_data.fdl_l = 0.43;
    sigma_data.flex_hal_l = 1.09;     sigma_data.fhl_l = 1.09;
    sigma_data.tib_ant_l = 1.91;      sigma_data.tibant_l = 1.91; 
    sigma_data.per_brev_l = 1.83/2;   sigma_data.perbrev_l = 1.83/2; % peroneals (brev/long) are combined
    sigma_data.per_long_l = 1.83/2;   sigma_data.perlong_l = 1.83/2; % peroneals (brev/long) are combined
    sigma_data.per_tert_l = 1.83/3;   % peroneals (tertius) assumed to be 1/3
    sigma_data.ext_dig_l = 0.43;     sigma_data.edl_l = 0.43;
    sigma_data.ext_hal_l = 1.09;     sigma_data.ehl_l = 1.09;
    sigma_data.ercspn_l = 99.0; %Back muscle not described
    sigma_data.intobl_l = 99.0; %Back muscle not described
    sigma_data.extobl_l = 99.0; %Back muscle not described    
    sigma_data.add_long_l = 2.26;    sigma_data.addlong_l = 2.26;
    sigma_data.add_brev_l = 1.47;    sigma_data.addbrev_l = 1.47;
    
    mus_vol = zeros(length(muscleNames),1);
    for i = 1:length(muscleNames)
        if isfield(sigma_data,muscleNames{i})
            mus_vol(i,1) = sigma_data.(muscleNames{i});
        else
            mus_vol(i,1) = NaN;
            disp(['Muscle volumen of ' muscleNames{i} ' not found in script']);
        end
    end
    
end   