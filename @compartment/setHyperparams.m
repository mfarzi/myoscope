function setHyperparams(obj, varargin)
    % setHyperparams is a method for class MULTICOMPARTMENT
    %
    %   setHyperparams(obj, p) update the class property "hyperparams" with the
    %   input vector p.
    %
    %   setHyperparams(obj, <paramName>, <paramValue>) update hyperparameters
    %   using params-value pairs. See getHyperparamsName for a full list of  
    %   hyperparameters.
    %
    
    % special case when no hyperparameters exist
    if obj.nHyperparams == 0
        assert(nargin==1 || (nargin==2 && isempty(varargin{1})), ...
            'MATLAB:compartment:setHyperparams', ...
            'Model %s has no hyperparameters to set.', obj.name);
        return;
    end
            
    % check the number of input arguments
    assert(nargin>1, 'MATLAB:compartment:setHyperparams', ...
        'Too few input parameters.');
    
    if nargin==2
        % setHyperparams(obj, p)
        p = varargin{1};
    else
        % setHyperparams(obj, <paramName>, <paramValue>)
        p = readPairs(obj, varargin{:});
    end

    validateattributes(p, {'numeric'},...
        {'column', 'nrows', obj.nHyperparams, 'integer', 'nonnegative'},...
        'setHyperparams', 'p');

    obj.hyperparams = p;
    
    if obj.nCompartments>1
        iStartHparams = 1;
        for n = 1:obj.nCompartments
            iEndHparams= iStartHparams+ obj.comp{n}.nHyperparams-1; 
            obj.comp{n}.setHyperparams(p(iStartHparams:iEndHparams));
            iStartHparams = iEndHparams+1;
        end
    end
end

function p = readPairs(obj, varargin)
    % helper function to parse input params-value pairs

    % default parameters
    [p, hyperparamsList] = obj.getHyperparams();

    % update specific parameters 
    for i=1:2:nargin-1
        % read hyperparams name
        thisParam = varargin{i};
        assert(isa(thisParam, 'char'),...
            'MATLAB:multicompartment:setHyperparams',...
            'Input hyperparameter name must be of type char');
        idx = strcmp(hyperparamsList, thisParam);
        assert(sum(idx)==1, 'MATLAB:compartment:setHyperparams',...
            strcat("'%s' is not a valid hyperparameter name.",...
            '\nValid hyperparameter names: %s.'),...
            thisParam, strjoin(hyperparamsList, ','));
        
        % read hyperparams value
        v = varargin{i+1};
        validateattributes(v, {'numeric'}, ...
            {'nonnegative', 'scalar', 'integer'},...
            'multicompartment.setHyperparams', thisParam);
        p(idx) = v;
    end
end