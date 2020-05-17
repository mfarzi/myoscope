function [scheme, gVec] = myoscopeReadScheme(path2scheme)
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
    gVec = [];
    nDir = 0;
    for n = 1:N
        thisVec = [scheme.x(n); scheme.y(n); scheme.z(n)];
        if sum(thisVec)==0
            continue;
        end
        
        thisVec = thisVec/norm(thisVec);
        if isempty(gVec)
            gVec = [gVec; thisVec'];
            nDir = 1;
            scheme.G_dir(n) = 1;
        else
            angleDiff = acos(gVec*thisVec)*180/pi;
            % The gradient directions could be flipped; -g and g are the 
            % same.
            idx = angleDiff>90;
            angleDiff(idx) = 180 - angleDiff(idx);
            [minAngleDiff, k] = min(angleDiff);
            if minAngleDiff>10
                gVec = [gVec; thisVec'];
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