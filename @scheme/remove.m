function remove(obj, idx)
    % remove is a public method for class SCHEME
    %
    % remove(obj, idx) remove measurments at given indices. 
    %
    %   input arguments:
    %               obj: A scheme object 
    %               idx: A boolean vector of size nMeasurements or a vector
    %                    of indices.
    %
    % see also: add, select
    %% return type: bvalue
    if strcmp(obj.type, 'bvalue')
        % update object 
        obj.ghat(idx) = [];
        obj.bval(idx) = [];
        
        if ~isempty(obj.ghatCode)
            obj.ghatCode(idx) = [];
        end
        
        if ~isempty(obj.bvalCode)
            obj.bvalCode(idx) = [];
        end
    end
    
    %% return type: stejskal-tanner
    if strcmp(obj.type, 'stejskal-tanner')
        % update object 
        obj.ghat(idx) = [];
        obj.gmag(idx) = [];
        obj.dtCode(idx) = [];
        obj.deltaCode(idx) = [];
        obj.teCode(idx) = [];
        obj.bval(idx) = [];
        
        if ~isempty(obj.ghatCode)
            obj.ghatCode(idx) = [];
        end
        
        if ~isempty(obj.bvalCode)
            obj.bvalCode(idx) = [];
        end
    end
end





