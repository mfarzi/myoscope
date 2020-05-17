function [directions, magnitudes, bvalMatrix] = procpar_to_gradients(path2procpar, varargin)
%procpar_to_gradient_directions returns gradient directions
    
    % parse the param-value pairs pass through varargin
    p = inputParser;
    p.CaseSensitive = false;
    errMsg = 'Value must be scalar and logical.';
    isNominalParamValid = @(v) assert(islogical(v) && isscalar(v), errMsg);
    p.addParameter('Nominal'   , false, isNominalParamValid);
    
    p.parse(varargin{:});
    returnNominalDirections = p.Results.Nominal;
    
    % read required parameters from procpar
    params = searchprocpar(path2procpar, 'bvalue', 'bvalrr', 'bvalpp', ...
                                          'bvalss', 'bvalrp', 'bvalrs', ...
                                          'bvalsp', 'gdiff' , 'dro'   , ...
                                          'dpe'   , 'dsl');
    
    nominalDirections = [params.dro;params.dpe;params.dsl]';
    scale = sqrt(sum(nominalDirections.^2, 2));
    magnitudes = params.gdiff/100*scale;
    
    if returnNominalDirections
        directions = nominalDirections./(scale+eps);
    else
        bvalLen = length(params.bvalrr);
        bvalMatrix=zeros(3,3,bvalLen);
        bvalMatrix(1,1,:) = params.bvalrr(:);
        bvalMatrix(2,2,:) = params.bvalpp(:);
        bvalMatrix(3,3,:) = params.bvalss(:);
        bvalMatrix(2,1,:) = params.bvalrp(:);
        bvalMatrix(1,2,:) = params.bvalrp(:);
        bvalMatrix(3,1,:) = params.bvalrs(:);
        bvalMatrix(1,3,:) = params.bvalrs(:);
        bvalMatrix(2,3,:) = params.bvalsp(:);
        bvalMatrix(3,2,:) = params.bvalsp(:);
        
        directions = zeros(bvalLen, 3);
        for i = 1:bvalLen
            [vec, val] = eig(bvalMatrix(:,:,i));
            [~, idx] = sort(diag(val), 'ascend');
            vec = vec(:,idx);
            directions(i, :) = vec(:,3);
        end
        
        % b-value is never precisely zero. Correct directions with very low
        % b-values.
        directions(magnitudes==0,:) = 0;
        %\\
    end % of if returnNominalDirections