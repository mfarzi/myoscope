function [scheme, G_dir] = camino_read_scheme(path2scheme)
fileID = fopen(path2scheme,'r');
keepReading = true;

while keepReading
    thisLine = fgetl(fileID);
    if isa(thisLine, 'integer') && thisLine == -1
        keepReading = false;
    end
    k = strfind(thisLine,'VERSION: ');
    if ~isempty(k)
        schemeType = thisLine(k+9:end);
        keepReading = false;
    end
end

if strcmp(schemeType, 'STEJSKALTANNER')
    scheme = fscanf(fileID, '%f');
    scheme = reshape(scheme, [7, length(scheme)/7]);
    scheme = scheme';
    scheme = array2table(scheme, 'VariableNames', ...
        {'x', 'y', 'z', 'G_mag', 'DELTA', 'delta', 'TE'});
    
    % add the direction
    N = size(scheme, 1);
    scheme.G_dir = zeros(N, 1);
    G_dir = [];
    nDir = 0;
    for n = 1:N
        this_dir = [scheme.x(n); scheme.y(n); scheme.z(n)];
        if sum(this_dir)==0
            continue;
        end
        
        if isempty(G_dir)
            G_dir = [G_dir; this_dir'];
            nDir = 1;
            scheme.G_dir(n) = 1;
        else
            dis = sum((G_dir - repmat(this_dir', nDir, 1)).^2, 2);
            k = find(dis<eps);
            if isempty(k)
                G_dir = [G_dir; this_dir'];
                nDir = nDir + 1;
                scheme.G_dir(n) = nDir;
            else
                scheme.G_dir(n) = k;
            end
        end
    end 
        
else
    scheme = [];
end