function [passiveParm]= loadPassiveParam(optParam,path)

    side_list={'r','l'}; % first r as in the model
    kpe_side=NaN(1,80);
    ksf_side=NaN(1,80);
    for side_opt=1:length(side_list)
        loadData=load(path);
        Results=loadData.sim.Results;
        
        kpe_side(side_opt,:) =Results.onlyParam.(optParam{:}).kpe;   ksf_side(side_opt,:) =Results.onlyParam.(optParam{:}).ksf;
    end
    
    % Assumed symmetry
    kpe=[kpe_side(1,1:40) kpe_side(1,1:40)]; % 40 muscle in a leg side of the model
    ksf=[ksf_side(1,1:40) ksf_side(1,1:40)];    
    ke0=ones(1,80)*0.6;
    
    passiveParm  = [kpe; ksf; ke0];
end