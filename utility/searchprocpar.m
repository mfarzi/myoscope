function params = searchprocpar(procpar_file,varargin)

% Search Varian procpar and extract as many variables as specified.
% 
% Example: 
% params=search_procpar(procpar_file,'np','nv','nv2');
% cellfun(@(n,v) assignin('caller',n,v),fieldnames(params),struct2cell(params));
% NB: The second line expands the 'params' structure into separate variables (Optional)
% 
% Created, I.Teh, 23 Feb 2007
% Convert output variables into structure, I.Teh, 27 Nov 2012


% Load procpar
fid = fopen(procpar_file);
procpar = textscan(fid, '%s');
fclose(fid);
procpar = procpar{1};

% Go to correct idx depending whether data is from Varian or not
tmp=strcmp(procpar,'decpat2');
[C I]=find(tmp);
if ~isempty(C)  %if Varian data
    idx_add=10;
else            % if not Varian data
    idx_add=0;
end
        
% Find parameters
for idx=1:length(varargin)
    variable=varargin{idx};
    tmp=strcmp(procpar,variable);
    [C I]=find(tmp);
    if ~isempty(C)
        C1=C(1)+1+idx_add;
        C2=C(1)+2+idx_add;
        num_elem=str2double(cell2mat(procpar(C1)));
        start=C2;
        for idx2=start:start+num_elem-1
            tmp=cell2mat(procpar(idx2));
            if strcmp(tmp(1),'"')
                value=tmp(2:length(tmp)-1);
            else
                value(idx2-start+1)=str2double(tmp);
            end
        end          
        varargout{idx}=value;
        clear value
    else
        varargout{idx}=[];
    end
end

% Convert strings to variables and assign values
for idx=1:length(varargin)
    eval(sprintf('params.%s=varargout{idx};',varargin{idx}));
end
