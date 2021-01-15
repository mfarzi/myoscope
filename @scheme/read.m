function [scheme, gVec] = read(path2scheme)
GAMMA = 2.6751525e8; % rad s-1 T-1
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
    
    nScheme = size(scheme, 1);
    
    % compute b-vallue 
    scheme.bval = (scheme.DELTA-scheme.delta/3).*((scheme.delta .*scheme.G_mag)*GAMMA).^2*1e-6;
    
    % compute nominal b-value
    classNo = zeros(nScheme, 1);
    bvalNominal = scheme.bval(1);
    classNo(1) = 1;
    nClass = 1;
    for n=2:nScheme
        tolerance = max(10, 0.1*bvalNominal);
        thisClassNo = find(abs(bvalNominal-scheme.bval(n))<tolerance);
        if isempty(thisClassNo)
            nClass = nClass + 1;
            classNo(n) = nClass;
            bvalNominal = [bvalNominal; scheme.bval(n)];
        else
            classSize = sum(classNo==thisClassNo);
            bvalNominal(thisClassNo) = (classSize*bvalNominal(thisClassNo)+...
                                        scheme.bval(n))/(classSize+1);
            classNo(n) = thisClassNo;
        end
    end
    scheme.bvalNominal = floor(bvalNominal(classNo));
    
    % add the direction
    N = size(scheme, 1);
    scheme.G_dir = zeros(N, 1);
    gVec = [];
    nDir = 0;
    for n = 1:N
        thisVec = [scheme.x(n); scheme.y(n); scheme.z(n)];
        if sum(thisVec)==0 || abs(scheme.bvalNominal(n))<10
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
    
    % consider repetitions
    % if DELTA, delta, TE, G_dir, and b-values are the same, then these
    % scans are considered as repeated measurements
    
    % consider different b-values with 10% variation
    tbl = [scheme.DELTA, scheme.delta, scheme.TE, scheme.G_dir, scheme.bvalNominal];
    [~,~,ic] = unique(tbl, 'rows');
    scheme.rep = ic;
else
    scheme = [];
end