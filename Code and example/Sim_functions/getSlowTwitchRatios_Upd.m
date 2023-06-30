% This function returns the percentage of slow twitch fibers in the muscles
% The data come from Uchida et al. (2016).
% We used 0.5 when no data were available.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% Updated: Israel Luis
% Date: 15/07/2022

function pctst = getSlowTwitchRatios_Upd(muscleNames)

    pctst_data.glut_med1_r = 0.55;      pctst_data.glmed1_r = 0.55;
    pctst_data.glut_med2_r = 0.55;      pctst_data.glmed2_r = 0.55;
    pctst_data.glut_med3_r = 0.55;      pctst_data.glmed3_r = 0.55;
    pctst_data.glut_min1_r = 0.55;      pctst_data.glmin1_r = 0.55;
    pctst_data.glut_min2_r = 0.55;      pctst_data.glmin2_r = 0.55;
    pctst_data.glut_min3_r = 0.55;      pctst_data.glmin3_r = 0.55;
    pctst_data.semimem_r = 0.4925;  % Yes yamaguchi
    pctst_data.semiten_r = 0.425;   % 50 yamaguchi
    pctst_data.bifemlh_r = 0.5425;      pctst_data.bflh_r = 0.5425;         pctst_data.bflh140_r = 0.5425; 
    pctst_data.bifemsh_r = 0.529;       pctst_data.bfsh_r = 0.529;          pctst_data.bfsh140_r = 0.529;  % 70 yamaguchi
    pctst_data.sar_r = 0.50;            pctst_data.sart_r = 0.50;
    pctst_data.add_mag1_r = 0.552;      pctst_data.addmagDist_r = 0.552;
    pctst_data.add_mag2_r = 0.552;      pctst_data.addmagIsch_r = 0.552;
    pctst_data.add_mag3_r = 0.552;      pctst_data.addmagMid_r = 0.552;     pctst_data.addmagProx_r = 0.552;
    pctst_data.tfl_r = 0.50;
    pctst_data.pect_r = 0.50;
    pctst_data.grac_r = 0.50;
    pctst_data.glut_max1_r = 0.55;      pctst_data.glmax1_r = 0.55;
    pctst_data.glut_max2_r = 0.55;      pctst_data.glmax2_r = 0.55;
    pctst_data.glut_max3_r = 0.55;      pctst_data.glmax3_r = 0.55;
    pctst_data.iliacus_r = 0.50;
    pctst_data.psoas_r = 0.50;
    pctst_data.quad_fem_r = 0.50;
    pctst_data.gem_r = 0.50;
    pctst_data.peri_r = 0.50;           pctst_data.piri_r = 0.50;
    pctst_data.rect_fem_r = 0.3865;     pctst_data.recfem_r = 0.3865; %0.32 Miller2014  0.45 yamaguchi
    pctst_data.vas_med_r = 0.503;       pctst_data.vasmed_r = 0.503;
    pctst_data.vas_int_r = 0.543;       pctst_data.vasint_r = 0.543;
    pctst_data.vas_lat_r = 0.455;       pctst_data.vaslat_r = 0.455;        pctst_data.vaslat140_r = 0.455; % 46.9 mean in 42 studies  - Human skeletal muscle fiber type percentage and area after reduced muscle use: A systematic review and meta-analysis
    pctst_data.med_gas_r = 0.566;       pctst_data.gasmed_r = 0.566;
    pctst_data.lat_gas_r = 0.507;       pctst_data.gaslat_r = 0.507;        pctst_data.gaslat140_r = 0.507;
    pctst_data.soleus_r = 0.803; % yes yamaguchi
    pctst_data.tib_post_r = 0.60;       pctst_data.tibpost_r = 0.60; % yes yamaguchi
    pctst_data.flex_dig_r = 0.60;       pctst_data.fdl_r = 0.60;  %50 yamaguchi
    pctst_data.flex_hal_r = 0.60;       pctst_data.fhl_r = 0.60;  %50 yamaguchi
    pctst_data.tib_ant_r = 0.70;        pctst_data.tibant_r = 0.70; % yes yamaguchi and Miller2014
    pctst_data.per_brev_r = 0.60;       pctst_data.perbrev_r = 0.60; % yes yamaguchi
    pctst_data.per_long_r = 0.60;       pctst_data.perlong_r = 0.60; % yes yamaguchi
    pctst_data.per_tert_r = 0.75;                                    % 35 yamaguchi
    pctst_data.ext_dig_r = 0.75;        pctst_data.edl_r = 0.75; %50 yamaguchi
    pctst_data.ext_hal_r = 0.75;        pctst_data.ehl_r = 0.75; %50 yamaguchi
    pctst_data.ercspn_r = 0.60;
    pctst_data.intobl_r = 0.56;
    pctst_data.extobl_r = 0.58;    
    pctst_data.add_long_r = 0.50;       pctst_data.addlong_r = 0.50;
    pctst_data.add_brev_r = 0.50;       pctst_data.addbrev_r = 0.50;
    
    
    pctst_data.glut_med1_l = 0.55;      pctst_data.glmed1_l = 0.55;      
    pctst_data.glut_med2_l = 0.55;      pctst_data.glmed2_l = 0.55;
    pctst_data.glut_med3_l = 0.55;      pctst_data.glmed3_l = 0.55;
    pctst_data.glut_min1_l = 0.55;      pctst_data.glmin1_l = 0.55;
    pctst_data.glut_min2_l = 0.55;      pctst_data.glmin2_l = 0.55;
    pctst_data.glut_min3_l = 0.55;      pctst_data.glmin3_l = 0.55;
    pctst_data.semimem_l = 0.4925;
    pctst_data.semiten_l = 0.425;
    pctst_data.bifemlh_l = 0.5425;      pctst_data.bflh_l = 0.5425;         pctst_data.bflh140_l = 0.5425;
    pctst_data.bifemsh_l = 0.529;       pctst_data.bfsh_l = 0.529;          pctst_data.bfsh140_l = 0.529;
    pctst_data.sar_l = 0.50;            pctst_data.sart_l = 0.50;
    pctst_data.add_mag1_l = 0.552;      pctst_data.addmagDist_l = 0.552;
    pctst_data.add_mag2_l = 0.552;      pctst_data.addmagIsch_l = 0.552;
    pctst_data.add_mag3_l = 0.552;      pctst_data.addmagMid_l = 0.552;     pctst_data.addmagProx_l = 0.552;
    pctst_data.tfl_l = 0.50;
    pctst_data.pect_l = 0.50;
    pctst_data.grac_l = 0.50;
    pctst_data.glut_max1_l = 0.55;      pctst_data.glmax1_l = 0.55;
    pctst_data.glut_max2_l = 0.55;      pctst_data.glmax2_l = 0.55;
    pctst_data.glut_max3_l = 0.55;      pctst_data.glmax3_l = 0.55;
    pctst_data.iliacus_l = 0.50;
    pctst_data.psoas_l = 0.50;
    pctst_data.quad_fem_l = 0.50;
    pctst_data.gem_l = 0.50;
    pctst_data.peri_l = 0.50;           pctst_data.piri_l = 0.50;
    pctst_data.rect_fem_l = 0.3865;     pctst_data.recfem_l = 0.3865;
    pctst_data.vas_med_l = 0.503;       pctst_data.vasmed_l = 0.503;
    pctst_data.vas_int_l = 0.543;       pctst_data.vasint_l = 0.543;
    pctst_data.vas_lat_l = 0.455;       pctst_data.vaslat_l = 0.455;        pctst_data.vaslat140_l = 0.455;
    pctst_data.med_gas_l = 0.566;       pctst_data.gasmed_l = 0.566;
    pctst_data.lat_gas_l = 0.507;       pctst_data.gaslat_l = 0.507;        pctst_data.gaslat140_l = 0.507;  
    pctst_data.soleus_l = 0.803;
    pctst_data.tib_post_l = 0.60;       pctst_data.tibpost_l = 0.60;
    pctst_data.flex_dig_l = 0.60;       pctst_data.fdl_l = 0.60;
    pctst_data.flex_hal_l = 0.60;       pctst_data.fhl_l = 0.60;
    pctst_data.tib_ant_l = 0.70;        pctst_data.tibant_l = 0.70; 
    pctst_data.per_brev_l = 0.60;       pctst_data.perbrev_l = 0.60;
    pctst_data.per_long_l = 0.60;       pctst_data.perlong_l = 0.60;
    pctst_data.per_tert_l = 0.75;
    pctst_data.ext_dig_l = 0.75;        pctst_data.edl_l = 0.75;
    pctst_data.ext_hal_l = 0.75;        pctst_data.ehl_l = 0.75;
    pctst_data.ercspn_l = 0.60;
    pctst_data.intobl_l = 0.56;
    pctst_data.extobl_l = 0.58;    
    pctst_data.add_long_l = 0.50;       pctst_data.addlong_l = 0.50;
    pctst_data.add_brev_l = 0.50;       pctst_data.addbrev_l = 0.50;
    
    pctst = zeros(length(muscleNames),1);
    for i = 1:length(muscleNames)
         if isfield(pctst_data,muscleNames{i})
            pctst(i,1) = pctst_data.(muscleNames{i});
         else
             pctst(i,1) = 0.5;
             disp(['Percentage of slow twitch fibers of ' muscleNames{i}...
                 '  not found in script getSlowTwitchRatios_lr. Used '...
                 'default value of 0.5']);
         end
    end    
end   