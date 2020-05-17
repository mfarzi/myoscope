function [image,k_space]=reconfid(path2fid, params, varargin)

% Reconstructs images from Varian .fid files
% Created as recon_cartesian.m, I.Teh, 31 Aug 2006
% Major revision to expand list of sequences reconstructed as recon_fid.m, I.Teh, 12 Feb 2014

% Amended by M Farzi 28 Nov 2018 
% _________________________________________________________________________
pList = {'removeDC', 'filtertype'};
     
nParameters = nargin - 2;
if nParameters > 0 
    for i = 1:2:nParameters
        iParameter = find(ismember(pList, varargin{i}));
        if isempty(iParameter)
            error(['%s is not a valid parameter.'                   ,...
                   '"%s", "%s", "%s", and "%s" are only allowed']   ,...
                   varargin{i}, pList{1} ,pList{2});
        end
        
        switch iParameter
            case 1
                remove_dc = varargin{i+1};
                if ~isa(remove_dc, 'logical')
                    error('Undefined parameter of type %s for %s\n',...
                           class(remove_dc), pList{1});
                end
            case 2
                filterType = varargin{i+1};
                if ~isa(filterType, 'char')
                    error('Undefined parameter of type %s for %s\n',...
                           class(filterType), pList{2});
                end
        end
    end
end

if ~exist('remove_dc', 'var')
    % Remove DC artifact
    remove_dc = true;
end

if ~exist('filterType', 'var')
    filterType = ''; % Available: '', 'butter'
end

% Setup directory and files to read
% if exist('directory','var')~=1 
%     disp ('Pick any .fid file in image dataset')
%     [~,PathName] = uigetfile;
%     directory=PathName(1:length(PathName)-7);
%     filestart = input('Specify first folder number (#) >> ');
%     fileend = input('Specify last folder number (#) >> ');
% else
%     filestart=file(1);
%     fileend=file(length(file));
% end

%for idx_file=filestart:fileend    
% Define write_dir
%str_k=num2str_zeros(idx_file,2);
%directory_full=strcat(directory,str_k);

% Initialise variables
cellfun(@(n,v) assignin('caller',n,v),fieldnames(params),struct2cell(params));     
ro=np/2;
pe=nv;
pe2=nv2;
if strcmp(seqfil_base,'epip')
    ro=nread/2;
    pe=nphase;
    pe2=nphase2;
end    
petable=lower(petable);
epi=strfind(seqfil,'epi');                             % EPI type sequence
epi3d=strfind(seqfil,'epi3d');                         % EPI3D type sequence
mems=strfind(seqfil,'mems');                            % ME-SE type sequence
mgems=strfind(seqfil,'mgems');                          % ME-GRE type sequence
if isempty(b0n), b0n=0; end
if isempty(b1n), b1n=0; end
if isempty(dro), dro=0; end
if isempty(dpe), dpe=0; end
if isempty(dsl), dsl=0; end
if nv2>0, ns=nv2; end
if isempty(etl), etl=1; end    
if isempty(params.diff), params.diff='n'; end         
if isempty(nphase2), nphase2=0; end       

% Load raw .fid data
[RE,IM,NP,NB,NT] = loadfid(path2fid);
k_space=complex(RE, IM);
clear RE; clear IM; 
nrcvrs=length(strfind(rcvrs,'y'));
image_n_all=size(k_space,3);

% 1D Butterworth filter (along readout)
if ~isempty(filterType)
    switch filterType
        case 'butter'
            f = filter_butter(ro,30,0.8)'; 
    end 
    siz=size(k_space);
    k_space = k_space.*repmat(f,[1 siz(2:end)]);               
end

% Reconstruct data
switch seqfil_base
    case {'mems','mgems'}    
        k_space=permute(k_space,[1,3,2]);                           
        image=zeros(size(k_space,1),size(k_space,2),size(k_space,3));
        for idx=1:size(k_space,3)
            image(:,:,idx)=fftshift(fft2(fftshift(k_space(:,:,idx))));
        end
        image=reshape(image,[size(image,1),size(image,2),params.ne,ns]);
        image=permute(image,[1,2,4,3]);

    case {'epi','epip'}  
        % epiref_type='manual';
        [image,k_space] = recon_EPIP(k_space,ro,pe,ns,nphase2,nseg,nnav,epiref_type,seqfil_type,params.orient,pescheme,params); 

    case 'se3d'
        k_space = reshape(k_space,[ro,pe,nv2]);
        image=fftshift(fftn(ifftshift(k_space)));

    case 'me3d'
        k_space = reshape(k_space,[ro,params.ne,pe,nv2]);
        k_space=permute(k_space,[1,3,4,2]);
        for idx=1:size(k_space,4)
            image(:,:,:,idx)=fftshift(fftn(ifftshift(k_space(:,:,:,idx))));
        end

    case 'fse3d' 
        if strcmp(cs_flag,'y')

            if length(idx_array)==1
                petable_array={petable_array};
                petable_array2={petable_array2};
            end

            for idx_rep=1:length(idx_array);
                idx_petable(:,:,idx_rep)=read_petable(petable_array{idx_rep})';
                idx_petable2(:,:,idx_rep)=read_petable(petable_array2{idx_rep})';
            end

            nvmin=min(idx_petable(:));
            nvmin2=min(idx_petable2(:));
            k=zeros(np/2,nvnom,nvnom2,length(idx_array));
            k_space=reshape(k_space,[np/2 nv nv2 length(idx_array)]);

            for idx_rep=1:length(idx_array)
                for idx_nv2=1:nv2
                    for idx_nv=1:nv
                        k(:,idx_petable(idx_nv,idx_nv2,idx_rep)-nvmin+1,idx_petable2(idx_nv,idx_nv2,idx_rep)-nvmin2+1,idx_rep)=k_space(:,idx_nv,idx_nv2,idx_rep);
                    end
                end
            end
            nv=nvnom;
            nv2=nvnom2;
            k_space_unsorted=k_space;
            k_space=k;
        else
            try
                k_space(:,pelist+pe/2,:)=k_space;
                k_space = reshape(k_space,[ro,pe,numel(dro),nv2]);
                k_space=permute(k_space,[1 2 4 3]);
            catch
                try
                    idx_petable=read_petable(petable);
                    if size(idx_petable,2)>etl %CS reconstruction
                        k_space_cs = zeros(np/2,nv,nv2);
                        % idx_petable = idx_petable(1:200,:); %tmp truncate
                        k_space = reshape(k_space,[np/2 size(k_space,2)*size(k_space,3) numel(k_space)/(np/2*size(k_space,2)*size(k_space,3))]);
                        idx_petable2 = repmat(idx_petable(:,1),[1 etl]);
                        idx_petable = idx_petable(:,2:end);
                        idx_petable = idx_petable';
                        idx_petable2 = idx_petable2';
                        idx_petable = idx_petable(:);
                        idx_petable2 = idx_petable2(:);
                        idx_petable = idx_petable-min(idx_petable)+1;
                        idx_petable2 = idx_petable2-min(idx_petable2)+1;
                        for idx=1:size(idx_petable,1)
                            k_space_cs(:,idx_petable(idx),idx_petable2(idx),:)=k_space(:,idx,:);
                        end
                        k_space=k_space_cs;
                    else
                        idx_petable=reshape(idx_petable',[1 numel(idx_petable)]);
                        k_space(:,idx_petable-min(idx_petable)+1,:,:)=k_space;
                    end
                catch
                    try                          
                        % IT: Ignore FSE3D navigator, 3 Jul 15
                        if strcmp(navigator,'y')
                            siz=[size(k_space) 1];
                            k_space_all = zeros([ro siz(2)/(etl+nnav)*etl siz(3)]);
                            phase_all = zeros([ro siz(2)/(etl+nnav)*nnav siz(3)]);
                            for idx_vol=1:siz(3)
                                k_space_curr=reshape(k_space(:,:,idx_vol),[ro etl+nnav siz(2)/(etl+nnav)]);
                                phase_curr=k_space_curr(:,etl+1:end,:);    
                                k_space_curr(:,etl+1:end,:)=[];    

                        %         % apply phase correction
                        %         S = fftshift(fft(fftshift(k_space_curr,1),size(k_space_curr,1),1),1);
                        %         N = fftshift(fft(fftshift(phase_curr,1),size(phase_curr,1),1),1);
                        %         N = repmat(N,[1 etl/nnav 1]);    
                        %         S = abs(S).*exp(i*(angle(S)-angle(N)));
                        %         k_space_curr = fftshift(ifft(fftshift(S,1),size(S,1),1),1);

                                phase_all(:,:,idx_vol)=reshape(phase_curr,[ro siz(2)/(etl+nnav)*nnav]);
                                k_space_all(:,:,idx_vol)=reshape(k_space_curr,[ro siz(2)/(etl+nnav)*etl]);
                            end
                            k_space = k_space_all;    
                            clear k_space_all
                        end

                        for idx_k=1:size(k_space,3)   
                            tmp=reshape(k_space(:,:,idx_k),[ro pe pe2]);
                            tmp2(:,pelist+pe/2,:,:)=tmp;
                            k_space2(:,:,:,idx_k)=tmp2;
                        end                        
                        k_space=k_space2;
                    catch
                        %for single echo
                        k_space=reshape(k_space,[ro pe pe2 numel(dro)]);
                    end
                end
            end
        end
        d4=max([numel(dro) numel(gain)]);
        image=single(zeros(ro,nv,nv2,d4));            
        for idx=1:d4
            image(:,:,:,idx) = fftshift(fftn(ifftshift(k_space(:,:,:,idx))));
        end         
        if strcmp(cs_flag,'y')
            k_space=k_space_unsorted;
        end                

    case 'ge2d'
        k_space(:,pelist+pe/2,:)=k_space; % there was a +1 here that Darryl removed 7 March 2017
        for idx=1:size(k_space,3)
            image(:,:,idx)=fftshift(fft2(fftshift(k_space(:,:,idx))));
        end

    case 'ge2dT1'
        ne0=params.ne;
        image=zeros(np/2,nv,ne0,rcvrs); %use ne0 because ne is Matlab reserved word
        for idx_rcvrs=1:size(k_space,3)
            k=reshape(k_space(:,:,idx_rcvrs),[np/2 etl ne0 nv/etl]);
            k=permute(k,[1 2 4 3]);
            k=reshape(k,[np/2 nv ne0]);
            for idx_ne=1:ne0
                k(:,pelist+pe/2+1,idx_ne)=k(:,:,idx_ne);
                image(:,:,idx_ne,idx_rcvrs)=fftshift(fft2(ifftshift(k(:,:,idx_ne))));
            end
        end

    case 'ge3d'
        idx_ge3d=reshape(1:size(k_space,3),numel(dro),size(k_space,3)/numel(dro));
        idx_ge3d=idx_ge3d';
        idx_ge3d=idx_ge3d(:);
        k_space=k_space(:,:,idx_ge3d);
        k_space=reshape(k_space,[ro pe nv2 numel(dro)]);
        image=single(zeros(size(k_space)));
        for idx=1:numel(dro)
            image(:,:,:,idx) = fftshift(fftn(ifftshift(k_space(:,:,:,idx))));
        end  
        image=image(:,end:-1:1,end:-1:1,:);


    case 'ge3dseg'
        k_space=reshape(k_space,[nro npe ns 2]);
        k_space(:,pelist-min(pelist)+1,:,:)=k_space;
        for idx=1:size(k_space,4)
            image(:,:,:,idx)=fftshift(fftn(ifftshift(k_space(:,:,:,idx))));
        end

    case 'epi3d'
        image=main_EPI_recon_3D(k_space,ro,pe,nv2,petable); 

    case 'radial'
        image=recon_2dnc(directory,file,params);

    case 'sems'            
        k_space = reshape(k_space,[nro ns numel(k_space)/nro/npe/ns npe]);
        k_space = permute(k_space,[1 4 2 3]);

        for idx_shot=1:size(k_space,4)
            for idx_slice=1:ns
                image(:,:,idx_slice,idx_shot)=fftshift(fft2(ifftshift(k_space(:,:,idx_slice,idx_shot))));
            end
        end

    case {'spuls','prescanpower'} %'spuls'  
        k_space=reshape(k_space,[size(k_space,1) size(k_space,2)*size(k_space,3)]);
        image=zeros(size(k_space));
        for idx_rep=1:size(k_space,2)
            image(:,idx_rep)=fftshift(fft2(ifftshift(k_space(:,idx_rep))));
        end
        image=reshape(image,[size(k_space,1)*size(k_space,2) 1]);

    otherwise
        if strcmp(seqfil,'gems'), k=k_space; end
        if strcmp(navigator,'y'), k_space=k_space(:,1:2:end,:); end %remove navigator

        if strcmp(propeller,'y')
            shots=numblades;
            tmp=(repmat(1:shots,[etl 1])-1)*etl;
            pelist_untouched = pelist;
            pelist=repmat(pelist,[1 shots])+tmp(:)';
            pe=etl;

            k=reshape(k_space,[ro etl ns shots image_n_all]);
            k=permute(k,[1 2 4 3 5]);
            k=reshape(k,[ro etl*shots ns image_n_all]);                        
            if etl>1 
                k(:,pelist+pe/2,:,:)=k;       
            end
            shots0=shots;

%             % Testing
%             k=reshape(k,[ro etl ns*shots*image_n_all]);
%             image=zeros(size(k));
%             for idx=1:size(k,3)
%                 image(:,:,idx)=fftshift(fft2(fftshift(k(:,:,idx))));
%             end

            % Exclude blades
            exclude=''; % exclude every nth shot
            if ~isempty(exclude)
                idx_exclude=1:exclude:shots;                        
                anglelist(idx_exclude)=[];
                shots=shots*(exclude-1)/exclude;
            end

            % Generate trajectory
            us=2*etl*numblades/pi/nv;
            [x y] = gen_samples(np/2,etl,[],'propeller',us,anglelist);





            % Process data
            image=zeros(ro,ro,ns,image_n_all);
            k_space=zeros(ro,ro,ns,image_n_all);
            for idx_image=1:image_n_all
                for idx_slice=1:ns
                    k_curr = k(:,:,idx_slice,idx_image);
                    k_curr = reshape(k_curr,[ro etl shots0]);

                    % Exclude blades
                    if ~isempty(exclude)
                        k_curr(:,:,idx_exclude)=[];          
                    end

                    % Phase correction (remove trianglular 2D imagephase)                        
                    win_ro = repmat(triang(ro),[1 etl]);
                    win_pe = repmat(triang(etl),[1 ro])';
                    kwin = win_ro.*win_pe;
                    kwin = kwin ./ max(kwin(:));

                    for idx=1:shots       
                        im_curr(:,:,idx)=fftshift(fft2(ifftshift(k_curr(:,:,idx))));
                        im_win(:,:,idx)=fftshift(fft2(ifftshift(k_curr(:,:,idx).*kwin)));     
                        im_curr(:,:,idx) = abs(im_curr(:,:,idx)).*exp(1i*(angle(im_curr(:,:,idx))-angle(im_win(:,:,idx))));
                        k_curr(:,:,idx)=fftshift(ifft2(ifftshift(im_curr(:,:,idx))));

                    end


                    % DMC 13/8/15 - create the T2 structure for T2w
                    % correction
                    T2_struct = {};

                    TE_mask = zeros(size(k_curr)); % DMC 13/8/15 - keep track of echoes
                    tmp = zeros(1,etl);
                    tmp(pelist_untouched + etl/2) = 1:etl;
                    for idx = 1:etl
                        TE_mask(:,idx,:) = tmp(idx);
                    end

                    TE = te + (0:(etl-1)) * esp;

                    TE_mask = reshape(TE_mask,[ro etl*shots]);

                    T2_struct.T2 = T2(:,:,idx_slice);
                    T2_struct.TE_mask = TE_mask(:);
                    T2_struct.TE = TE;


                    % Regrid
                    k_curr = reshape(k_curr,[ro etl*shots]);


                    [image(:,:,idx_slice,idx_image),k_space(:,:,idx_slice,idx_image)] = nufft_propeller(x,y,k_curr(:)',np/2,nv, T2_struct);

%                        [image(:,:,idx_slice,idx_image),k_space(:,:,idx_slice,idx_image)] = nufft_propeller(x,y,k_curr(:)',np/2,nv); %NUFFT
%                         [image(:,:,idx_slice,idx_image),k_space(:,:,idx_slice,idx_image)] = kb_2d(x,y,k_curr(:)',ro); %Kaiser-Bessel
                end
            end
        else
            k=reshape(k_space,[ro etl ns shots image_n_all]);
            k=permute(k,[1 2 4 3 5]);
            k=reshape(k,[ro etl*shots ns image_n_all]);   
            if ~isempty(pelist)
                k(:,pelist+pe/2,:,:)=k;  
            end

            % 2D FFT
            image=single(zeros(size(k)));
            for image_n=1:size(k,4)
                for slice=1:size(k,3)
                    image(:,:,slice,image_n)=fftshift(fft2(fftshift(k(:,:,slice,image_n))));
                end
            end     
        end

        % Sum of squares addition for multiple receive channels
        if nrcvrs>1
            image_n_all=image_n_all/nrcvrs;
            for slice=1:ns       
                for image_n=1:image_n_all
                    idx=nrcvrs*(image_n-1)+1:nrcvrs*image_n;
                    image(:,:,slice,image_n)=sum(abs(image(:,:,slice,idx)).^2,4).^0.5;                        
                end                        
            end
            image(:,:,:,image_n_all+1:size(image,4))=[];
        end
end    

% Reorder according to slice order (eg. linear or interleaved) if 2D multislice
try
    if (isempty(nv2) || nv2==0) && ~strcmp(params.orient,'3orthogonal') && ns>1
        [~,idx]=sort(pss);
        image=image(:,:,idx,:);  
    end 

    % Remove DC artifact
    if remove_dc
        if nv2>0 || nphase2>0 || nvnom2>0
            if strcmp(cs_flag,'y')
                pe=nvnom;
                pe2=nvnom2;
            end

            %Find peak index in centre 3x3x3 neighbourhood
            ctr=image(ro/2-1:ro/2+1,pe/2-1:pe/2+1,pe2/2-1:pe2/2+1,1);
            [C tmp]=max(ctr);
            [C2 tmp]=max(squeeze(C));
            [tmp I3]=max(squeeze(C2));    
            [tmp I2]=max(squeeze(C(:,:,I3)));
            [tmp I1]=max(squeeze(ctr(:,I2,I3)));
            I1=ro/2-2+I1;
            I2=pe/2-2+I2;
            I3=pe2/2-2+I3;            

            for idx=1:size(image,4)
                %Replace peak with mean of voxels surrounding peak
                voxs=1;
                tmp=image(I1-voxs:I1+voxs,I2-voxs:I2+voxs,I3-voxs:I3+voxs,idx);
                image(I1,I2,I3,idx)=mean(tmp(:));                    
                k_space(:,:,:,idx) = fftshift(ifftn(ifftshift(image(:,:,:,idx))));
            end
        end  
    end
end
%end
end

% *******************************************************************

function [image,k_out]=recon_EPIP(k,ro,pe,ns,nphase2,nseg,nnav,epiref_type,seqfil_type,orient,pescheme,params) 

% Phase correction preprocessing for EPIP data
% Adapted from Main_epi_recon_2D by B.Zheng

% User entered parameters
% tkw = 0.1; % Tukey window filter in nv2 direction (0=boxcar to 1=Hanning) Set to 0.1 for s_20130813_02

% Initialise parameters
switch epiref_type
    case 'manual'
        %[1,1,1,1,...]
        i0=[];              % i0 - Normal readout, no phase encoding
        iminus2=[];         % i-2 - Inverted readout, no phase encoding
        iminus1=[];         % i-1 - Inverted readout, with phase encoding
        i1=1:size(k,3);   % i1 - Normal readout, with phase encoding
    case 'single'
        %[0,1,1,1,1,...]
        i0=1;
        iminus2=[];
        iminus1=[];
        i1=(2:size(k,3));
    case 'triple'
        %[0,-2,-1,1,1,1,1,...]
        i0=1;
        iminus2=2;
        iminus1=3;
        i1=(4:size(k,3));
    case 'fulltriple'
        %[0,-2,-1,1,-1,1,-1,1,-1,1,...]
        i0=1;
        iminus2=2;
        iminus1=3:2:size(k,3);
        i1=4:2:size(k,3);
end
image_n_all=numel([i0 iminus2 iminus1 i1]);
etl=pe/nseg;

%Check for readout oversampling
if size(k,1)>=2*ro*etl
    over=2;
else
    over=1;
end
    
% Extract navigator data - Do something with navigator? 
if ~isempty(nnav)
    if nnav>0
        nav=k(1:(1+nnav)*over*ro,:,:);
        k(1:(1+nnav)*over*ro,:,:)=[];
    else
        nav=k(1:(1)*over*ro,:,:);
        k(1:(1)*over*ro,:,:)=[];
    end
end

% Order by readout x phase encode x slices x shots x images
if strcmp(seqfil_type,'3d')
    ns=nphase2;
end
k=reshape(k,[ro*over etl ns nseg image_n_all]);

% % Consider Partial Fourier
% pf=0.6; % Acquired k-space fraction
% pf=round((1-pf)*nphase2);
% k(:,:,end-pf+1:end,:,:)=conj(k(:,:,1:pf,:,:));

% Apply Tukey window (get rid of asymmetric k-space spike in nv2) 
if exist('tkw','var')
    w = tukeywin(nphase2,tkw); 
    w = repmat(w, [1 size(k,1) size(k,2) size(k,4) size(k,5)]);
    w = permute(w,[2 3 1 4 5]);
    k = k.*w;
    clear w
end

% FFT along slice for 3D
if strcmp(seqfil_type,'3d')
    k = fftshift(fft(fftshift(k,3),size(k,3),3),3);
end

% Reverse every second echo
k(:,2:2:end,:,:,:) = k(end:-1:1,2:2:end,:,:,:);  

% Reverse inverted polarity data
k(:,:,:,:,[iminus2 iminus1]) = k(end:-1:1,:,:,:,[iminus2 iminus1]);  

% Apply 1D fft in readout
s = fftshift(fft(fftshift(k,1),size(k,1),1),1);

if (strcmp(params.rf1off,'y') && strcmp(params.rf2off,'y'))||(params.flip1==0 && (params.flip2==0)) % noise measurement
    s = s(:,:,:,:,1);
    i1 = 1;
else
    % epiref_type='single';
    switch epiref_type
        case 'manual'
            % No phase correction. Use for noise measurements.
        case 'single'
            % Phase correction using non-phase encoded data
            s(:,:,:,:,i1)=abs(s(:,:,:,:,i1)).*exp(-i*(angle(s(:,:,:,:,i1))-repmat(angle(s(:,:,:,:,i0)),[1 1 1 1 numel(i1)])));
        case 'triple' % Needs work
            % Phase correction using non-phase encoded then inverted readout data
            s(:,:,:,:,i1)=abs(s(:,:,:,:,i1)).*exp(-i*(angle(s(:,:,:,:,i1))-repmat(angle(s(:,:,:,:,i0)),[1 1 1 1 numel(i1)])));
            s(:,:,:,:,iminus1)=abs(s(:,:,:,:,iminus1)).*exp(-i*(angle(s(:,:,:,:,iminus1))-angle(s(:,:,:,:,iminus2))));            
            phasediff=angle(s(:,:,:,:,i1)./repmat(s(:,:,:,:,iminus1),[1 1 1 1 numel(i1)]));            
            s(:,:,:,:,i1)=abs(s(:,:,:,:,i1)).*exp(-i*(angle(s(:,:,:,:,i1))-phasediff/2));
        case 'fulltriple'
            % Phase correction using non-phase encoded then inverted readout data
            s(:,:,:,:,i1)=abs(s(:,:,:,:,i1)).*exp(-i*(angle(s(:,:,:,:,i1))-repmat(angle(s(:,:,:,:,i0)),[1 1 1 1 numel(i1)])));
            s(:,:,:,:,iminus1)=abs(s(:,:,:,:,iminus1)).*exp(-i*(angle(s(:,:,:,:,iminus1))-repmat(angle(s(:,:,:,:,iminus2)),[1 1 1 1 numel(i1)])));
            phasediff=angle(s(:,:,:,:,i1)./s(:,:,:,:,iminus1));
            s(:,:,:,:,i1)=abs(s(:,:,:,:,i1)).*exp(-i*(angle(s(:,:,:,:,i1))-phasediff/2));
            s(:,:,:,:,iminus1)=abs(s(:,:,:,:,iminus1)).*exp(-i*(angle(s(:,:,:,:,iminus1))+phasediff/2));
            s(:,:,:,:,i1)=(s(:,:,:,:,i1)+s(:,:,:,:,iminus1))/2; % Average signals
        otherwise
           disp('Only single and full triple mode supported');
    end
end

% Trim readout border
if over==2
    s=s(ro/2+1:ro*3/2,:,:,:,:);
end

% Apply 1D ifft in readout
k = fftshift(ifft(fftshift(s,1),size(s,1),1),1);
clear s

% Reorder multishot data
if strcmp(pescheme,'c') && nseg>1
    pe_steps=0:nseg/2:pe/2-1;
else
    pe_steps=0:nseg:pe-1;
end
k_out=zeros(size(k,1),size(k,2)*size(k,4),size(k,3),numel(i1));
for idx_nseg=1:nseg
    if strcmp(pescheme,'c') && nseg>1
        idx1=sign(rem(idx_nseg,2)-0.5);
        idx0=pe/2+1+ceil((idx_nseg-1)/2)*idx1;
        k_out(:,idx1*pe_steps+idx0,:,:)=k(:,:,:,idx_nseg,i1);
    else
        k_out(:,pe_steps+idx_nseg,:,:)=k(:,:,:,idx_nseg,i1); 
    end    
end
clear k

% 2D FFT
image=single(zeros(size(k_out)));
for image_n=1:size(k_out,4)
    for slice=1:size(k_out,3)
        image(:,:,slice,image_n)=fftshift(fft2(fftshift(k_out(:,:,slice,image_n))));
    end
end

end

% *******************************************************************

function [image varargout]=recon_2dnc(path_full,file,params)

% Reconstructs 2D non-cartesian images from fid data
% 
% Created by I.Teh, 22 Jul 2013
% NEED TO FIX GRADIENT TIMING CORRECTION

%__________________________________________________________________________


% Setup directory and files to read
path_full=strcat(path_full,num2str_zeros(file,2),'.fid');
             
% Initialise parameters
cellfun(@(n,v) assignin('caller',n,v),fieldnames(params),struct2cell(params));
if ~exist('bidirectional','var'), bidirectional='n'; end

% Read data
[Re Im ro num_blocks] = load_fid(path_full(1:end-4));
k=complex(Re, Im);
clear Re, clear Im;

% Generate sampling
nv2=1; %2D for now
acq = 'radial';
us = 1; % Undersampling factor 
[x y z] = gen_samples(np/2,nv,nv2,acq,us);

if strcmp(bidirectional,'y')
    x=[x x];
    y=[y y];
    z=[z z];
    nd=2;
else
    nd=1;
end
shots=nv/etl*nd;
    
% Rearrange sampling based on pelist
x=x(:,pelist+1);
y=y(:,pelist+1);
x=x(:);
y=y(:);

% % Plot sampling (Diagnostic only)
% close all
% figure
% tmp=max([abs(x); abs(y)]);
% xlim([-tmp tmp]);
% ylim([-tmp tmp]);
% axis image
% hold on
% for idx=1:nv*nd
%     plot(x((idx-1)*nv+1:idx*nv),y((idx-1)*nv+1:idx*nv),'.');
%     pause(0.01);
% end

% 2D reconstruction
image=zeros(np/2,nv,ns);
for idx_ns=1:ns
    %echo-shot-slice
    % k_curr=k(:,(idx_ns-1)*nv+1:idx_ns*nv); 
    %echo-slice-shot    
    tmp=(repmat((idx_ns-1)*etl+(1:etl),[shots 1]) + repmat(0:etl*ns:size(k,2)-1,[etl 1])')'; 
    tmp=tmp(:); 
    k_curr=k(:,tmp); 
    
    % Correct for gradient timing [T. Block, ISMRM, 2011]
    if strcmp(bidirectional,'y')
        k_curr1=abs(k_curr(:,1:nv));
        k_curr2=flipud(abs(k_curr(:,nv+1:end)));
        s_curr1=fftshift(fft(ifftshift(k_curr1,1),[],1),1);
        s_curr2=fftshift(fft(ifftshift(k_curr2,1),[],1),1);
        % Find cross correlation function
        g=s_curr1.*conj(s_curr2);
        cc=fftshift(ifft(ifftshift(g,1),[],1),1);
        
        % Find region of support 
        [C I]=max(abs(k_curr1(:,1)));
        [C I]=find(abs(k_curr1(:,1))>C*0.1);        
        
        % x_shift based on 0 and 180 degree spokes (need to generalise) 
        tmp = angle(g(C,1));
        p = polyfit(1:length(tmp),tmp',1);
        x_shift=-p(1)*nv/(2*pi);
        
        % Find y_shift based on 90 and 270 degree spokes (need to generalise) 
        tmp = angle(g(C,nv/2+1));
        p = polyfit(1:length(tmp),tmp',1);
        y_shift=-p(1)*nv/(2*pi);
        
        % Find and apply shift for each spoke
        angle_step=pi/nv;
        phi_array=(pelist(1:nv)+1)*angle_step;
        delta_k=zeros(1,nv);
        for idx_nv=1:nv
            phi=phi_array(idx_nv);
            delta_k(idx_nv)=((cos(2*phi)+1)*x_shift+(-cos(2*phi)+1)*y_shift)/2;
        end
        
        % Apply correction
        slope=linspace(-pi,pi,np/2+1)'; 
        slope(1)=[];
        s_curr1=fftshift(fft(ifftshift(k_curr(:,1:nv),1),[],1),1);
        for idx_nv=1:nv
            s_curr1(:,idx_nv) = abs(s_curr1(:,idx_nv)).*exp((angle(s_curr1(:,idx_nv))+delta_k(idx_nv)*slope)*j);            
        end
        k_curr=fftshift(ifft(ifftshift(s_curr1,1),[],1),1);
    end
    
    k_curr=k_curr(:,1:nv);
    k_curr=flipud(k_curr);
    x=x(1:numel(k_curr));
    y=y(1:numel(k_curr));
%     image(:,:,idx_ns) = nufft_propeller(y,x,k_curr(:)',np/2,nv); %NUFFT
    image(:,:,idx_ns) = kb_2d(x,y,k_curr(:),ro); %Kaiser-Bessel
end
close all
display_multigeneric(abs(image));
axis image
colormap gray

end

% *******************************************************************

function [image,k_space] = kb_2d(x,y,k,ro);

scale=max([max(abs(x(:))) max(abs(y(:)))])*2;
x0=x./scale./(1+eps);
y0=y./scale./(1+eps);
locs=complex(x0,y0);
dcf = calcdcflut(locs,ro);
k_space = gridkb(locs,k,dcf,ro,2.5,2.5); 
image = fftshift(fft2(fftshift(k_space))); 

end

% *******************************************************************

function f = filter_butter(pts,order,cutoff);

[z,p,k]=buttap(order); 
[num,den]=zp2tf(z,p,k); 
wc=pi*cutoff;
[num,den]= lp2lp(num,den,wc); 
w=linspace(0,pi,pts/2);
H=freqs(num,den,w); 
f=abs(H);
f=[fliplr(f) f];

end

% *******************************************************************



% Modification History
% 23 Feb 2007, I. Teh   Added display images normalised to specified noise
% 23 Feb 2007, I. Teh   Cleaned up for SNR and DWI analysis
% 22 Jun 2007, I. Teh   Fixed multislice reconstruction code
% 19 Jul 2007, I. Teh   Fixed multiplane code
% 11 Sep 2007, I. Teh   Added 3D recon
% 26 Sep 2007, I. Teh   Debugged 3D recon 
% 14 Nov 2007, I. Teh   Added in catenate slices / images for display
% 9 Jan 2008,  I. Teh   Added EPI reconstruction
% 18 Jan 2008, I. Teh   Added large file splitting function. See 'load_fid_by_parts.m'
% 25 Mar 2008, I. Teh   Added plot mean signal variation vs DW directions
% 16 Apr 2008, I. Teh   Added mems recon for T2 mapping
% 25 Jul 2013, I. Teh   Added recon_EPIP, radial FSE recon and tidied up
% 12 Feb 2014, I. Teh   Added recon for ge3dseg, fse3d, mems, me3d, spuls
% 4 Aug 2015, I. Teh    Added PROPELLER
% 13 Aug 2015, D. McC   Broke PROPELLER
% 28 Nov 2018, M. Farzi Modified the input arguments