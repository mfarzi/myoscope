function params=readprocpar(path2procpar)

% Read procpar to setup reconstruction 
% Add new variables to read in Items A & B
% Created I.Teh 21 Nov 2012

% Amended by M Farzi 28 Nov 2018

    params=searchprocpar(path2procpar,...
        'seqfil','petable','pelist','pe_order','sampling','navigator','array','rcvrs','gain','ident',...                    
        'field_strength','orient','psf','v','v_str','ssc','simulation','seedarr','randseed',...                              	
        'nseg','ns','np','nv','nv2','nv3','lro','lpe','lpe2','thk','te','etl','esp','pss','effTE','nsa','ne','nt',...       
        'dwi','b0n','b1n','diffscheme','dro','dpe','dsl','diffpol','gdiff','tdelta','tDELTA',...                                       
        'propeller','petable_split','philips_dwrot','nblades','undersample_pe3','undersample_blades',...    
        'epiref_type','nread','nphase','nphase2','nnav','ss','zpe_flag',...                                                  
        'bidirectional','diff','psi','phi','theta','pescheme',...
        'rf1off','rf2off','rf3off','flip1','flip2','flip3','nshells',...
        'cs_flag','petable_array','petable_array2','nvnom','nvnom2','idx_array',...                                             
        'bvalrr','bvalpp','bvalss','bvalrp','bvalrs','bvalsp','bvalue','cs_speedup','max_bval','mean_bval','image',...  
        'pro','ppe','ppe2','nnav','dro2','dpe2','dsl2','anglelist','numblades','bladewiden'...
        );   
    
    params.nro=params.np/2;
    if isempty(params.etl)
        params.etl=1; 
    end
    params.npe=params.nv;
    params.shots=params.nv/params.etl;                % # shots
    params.rcvrs=length(strfind(params.rcvrs,'y'));   % # receivers
    if isempty(params.nnav)
        params.nnav=0;
    end
    
    % Find length of arrays
    params.image_n_all=1;
    try
        if ~isempty(params.array)
            params.array(~isletter(params.array))=' ';
            params.array=regexp(params.array, '([^ ]*)', 'match');
            tmp=search_procpar(path2procpar,params.array{1});
            tmp_name=fieldnames(tmp);
            params.image_n_all=length(getfield(tmp,tmp_name{1}));
        end
    end

    % Initialise base pulse sequence
    k=strfind(params.seqfil,'_');
    tmp=params.seqfil;
    tmp(k)=[];
    seqfil_base_array={'sems','semsdw','se3d','fsems','fse3d','epi','epip','me3d',...
        'radial','gems','ge2d','ge2dT1','ge3d','ge3dseg','mems','mgems','spuls','simulation','prescanpower'};
    for idx=length(seqfil_base_array):-1:1
        k=strfind(tmp,seqfil_base_array{idx});
        if ~isempty(k)
            params.seqfil_base=seqfil_base_array{idx};
            break
        end
    end 
    
    % Initialise sequence type based on seqfil (ie. 2D or 3D)
    k=strfind(params.seqfil,'3d');
    if ~isempty(k)
        params.seqfil_type='3d';
    else
        params.seqfil_type='2d';
    end
    
    % EPI only
    if ~isempty(params.nread)
        params.np=params.nread;
        params.nv=params.nphase;
        params.nv2=params.nphase2;
    end
    
    % 1D
    if isempty(params.nv)
        params.nv=1;
    end

    % Get dimensions
    if strcmp(params.cs_flag,'y')
        params.x_res = params.lro*10/(params.np/2);
        params.y_res = params.lpe*10/params.nvnom;
        params.x_dim = params.lro ./ params.x_res * 10;
        params.y_dim = params.lpe ./ params.y_res * 10;
        if params.nv2>0 %3D
            params.z_res = params.lpe2*10/params.nvnom2;
            params.z_dim = params.lpe2 ./ params.z_res * 10;
        else %2D
            params.z_res = params.thk; 
            params.z_dim = params.ns;
        end 
    else
        params.x_res = params.lro*10/(params.np/2);
        params.y_res = params.lpe*10/params.nv;
        params.x_dim = params.lro ./ params.x_res * 10;
        params.y_dim = params.lpe ./ params.y_res * 10;
        if params.nv2>0 %3D
            params.z_res = params.lpe2*10/params.nv2;
            params.z_dim = params.lpe2 ./ params.z_res * 10;
        else %2D
            params.z_res = params.thk; 
            params.z_dim = params.ns;
        end 
    end
    
    % 1D
    if isempty(params.bladewiden)
        params.bladewiden=1;
    end
    
    
end