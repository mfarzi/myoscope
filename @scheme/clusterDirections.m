function [ghatDic, ghatCode] = clusterDirections(ghat, TOL)
    % clusterDirections is a static method for class scheme
    
    if nargin==1
        TOL = 10;
    end
    
    assert(isscalar(TOL)&&isnumeric(TOL)&&TOL>0&&TOL<90,...
        'MATLAB:scheme:invalidInputArgument',...
        'Tolerance must be a positive scalar number between 0 and 90.');
    
    ghat = scheme.normaliseDirections(ghat);
    
    % online k-means clustering
    nMeasurements = size(ghat, 1);
    ghatCode = ones(nMeasurements, 1);
    ghatDic = [0 0 0];
    clusterSize = 0;
    nWords = 1;
    %
    for i = 1:nMeasurements
        this_g = ghat(i, :);
        if all(this_g==0)
            clusterSize(1) = clusterSize(1) + 1;
            continue;
        end
        
        theta = acosd(abs(this_g*ghatDic'));
        [minVal, iMin] = min(theta);
        if minVal<TOL
            % update the cluster centre
            N = clusterSize(iMin);
            this_dicWord = ghatDic(iMin,:);
            % rotate the gradient vector if required
            if (this_g*this_dicWord')<0
                this_g = -this_g;
            end
            ghatDic(iMin,:) = (this_dicWord*N+this_g)/(N+1);
            clusterSize(iMin) = N+1;
            ghatCode(i) = iMin;
        else
            % add a new cluster
            ghatDic = [ghatDic; this_g];
            clusterSize = [clusterSize; 1];
            nWords = nWords+1;
            ghatCode(i) = nWords;
        end
    end
    
    if clusterSize(1)==0
        % remove direction [0 0 0]
        ghatDic(1,:) = [];
        ghatCode = ghatCode-1;
    end
end     