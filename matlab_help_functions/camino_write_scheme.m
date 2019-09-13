function camino_write_scheme(thisScheme, path2scheme, varargin)
pList = {'scheme', 'ns', 'np', 'nf', 'nd'};
nParameters = nargin - 2;
if nParameters > 0 
    for i = 1:2:nParameters
        iParameter = find(ismember(pList, varargin{i}));
        if isempty(iParameter)
            error(['%s is not a valid parameter.'                   ,...
                   '"%s", "%s", "%s", and "%s" are only allowed']   ,...
                   varargin{i}, pList{1} ,pList{2}, pList{3},        ...
                   pList{4}, pList{5});
        end
        
        switch iParameter
            case 1
                schemeType = varargin{i+1};
                if ~isa(schemeType, 'char')
                    error('Undefined parameter of type %s for %s\n',...
                           class(schemeType), pList{1});
                end
                
                if strcmp(schemeType, 'STEJSKALTANNER') || strcmp(schemeType, 'BVECTOR')
                    error(['Unknown scheme type %s. Either "BVECTOR" ',...
                           'or "STEJSKALTANNER" is only supported.\n'],...
                           schemeType);
                end

            case 2
                ns = varargin{i+1};
                if ~isa(ns, 'double')
                    error('Undefined parameter of type %s for %s\n',...
                           class(ns), pList{2});
                end
                
            case 3
                np = varargin{i+1};
                if ~isa(np, 'double')
                    error('Undefined parameter of type %s for %s\n',...
                           class(np), pList{3});
                end
                
            case 4
                nf = varargin{i+1};
                if ~isa(nf, 'double')
                    error('Undefined parameter of type %s for %s\n',...
                           class(nf), pList{4});
                end
                
            case 5
                nd = varargin{i+1};
                if ~isa(nd, 'double')
                    error('Undefined parameter of type %s for %s\n',...
                           class(nd), pList{5});
                end

            otherwise
                error(['%s is not a valid parameter.'              ,...
                   '"%s", "%s", "%s", and "%s" are only allowed']  ,...
                   varargin{i}, pList{1} ,pList{2}, pList{3},       ...
                   pList{4}, pList{5});
        end
    end
end

if ~exist('schemeType', 'var')
    if size(thisScheme, 2)==7
        schemeType = 'STEJSKALTANNER';
    elseif size(thisScheme, 2)==3
        schemeType = 'BVECTOR';
    else
        error(['The type of the scheme file cannot be automatically ',...
               'identified.\n']);
    end
end

note = '#image dimension: ';
if exist('ns', 'var')
    note = strcat(note, sprintf('ns=%d ,', ns));
else
    note = strcat(note, 'ns=NA ,');
end

if exist('np', 'var')
    note = strcat(note, sprintf('np=%d ,', np));
else
    note = strcat(note, 'np=NA ,');
end
 
if exist('nf', 'var')
    note = strcat(note, sprintf('nf=%d ,', nf));
else
    note = strcat(note, 'nf=NA ,');
end

if exist('nd', 'var')
    note = strcat(note, sprintf('nd=%d', nd));
else
    note = strcat(note, 'nd=NA');
end

fileID = fopen(path2scheme,'w');
fprintf(fileID, sprintf('#%s scheme.\n', schemeType));
fprintf(fileID, strcat(note, '\n'));

if strcmp(schemeType, 'STEJSKALTANNER')
    fprintf(fileID, ['#x         y         z         |G|       DELTA     ',...
                     'delta     TE       \n']);
    fprintf(fileID, 'VERSION: STEJSKALTANNER\n');
    nScheme = size(thisScheme, 1);

    for i=1:nScheme
        A = thisScheme{i,1:7};
        for j=1:7
            if j<7 && A(j)<0
                fprintf(fileID, '%1.5f  ',A(j));
            elseif j<7 && A(j)>=0
                fprintf(fileID, ' %1.5f  ', A(j));
            elseif j==7 && i<nScheme
                fprintf(fileID, ' %1.5f\n', A(j));
            else
                fprintf(fileID, ' %1.5f', A(j));
            end
        end
    end
end

fclose(fileID);