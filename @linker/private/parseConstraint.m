function config = parseConstraint(obj, str)
    % parseConstraint is a private method for class LINKER.
    % 
    % parseConstraint(obj, str) parse the input string into a struct for
    % use by addConstraint.
    %       Input Arguments:
    %                   obj: Multicompartment object
    %                   str: Input of type char. Semantic parameter names
    %                        must be used. See getParamsName(obj) for a 
    %                        list of valid semantic names.
    %                        Examples:
    %                           1) 's0 = 1' set s0 as constant with value 1
    %                           2) 'diffPar >= 1.5e-9'
    %                           3) 'diffPar <= 1.5e-9'
    %                           4) 'diffPar >= diffPerp1'
    %                           5) 'diffPerp2 = diffPerp1'
    %
    %      Output Arguments:
    %                config: A struct with the following fields
    %                         str: input constraint
    %                        type: 'eq', 'ge', 'le'
    %                         lhv: A scalar positive integer number to 
    %                              indicate the left-hand variable
    %                         rhV: A positive integer vector of  the 
    %                              right-hand variables
    %                          fh: function handle that relates right-hand
    %                              variables to the left-hand ones
    %                          gh: function handle to the gradient of
    %                              right-hand expression wrt p
    
    % check input object
    validateattributes(obj, {'linker'}, {'column'});
    objLen = size(obj, 1);
    
    % check if input is of type char
    assert(ischar(str), 'MATLAB:LINKER:parseConstraint',...
        "Input constraint must be of type 'char'.");
    
    % remove all spaces
    str(isspace(str))=[];
    
    % indentify the left-hand variable
    varNames = obj.getName();
    lhv = cellfun(@(c) startsWith(str, c), varNames);
    assert(sum(lhv)==1, 'MATLAB:LINKER:parseConstraint', ...
        strcat('Input constraint must start with a valid semantic',...
        ' name including:\n%s.'), strjoin(varNames, ', '));
    config.lhv = find(lhv);    
    
    % trim str to exclude variable name
    varLen = length(varNames{lhv});
    str = str(varLen+1:end);
    
    % identify the type of constraint
    constraintList = {'=', '>=', '<='};
    type = cellfun(@(c) startsWith(str, c), constraintList);
    assert(sum(type)==1, 'MATLAB:LINKER:parseConstraint', ...
        strcat('Input constraint must contain a valid operator',...
        ' including: %s.'), strjoin(constraintList, ', '));
    config.type = constraintList{type};
    
    % trim str to exclude the operator
    varLen = length(constraintList{type});
    str = str(varLen+1:end);
    
    % assess the right-hand expression 
    % option 1: it is a number
%     [val, isaNumber] = str2num(str);
%     if isaNumber
%         config.rhv = [];
%         config.fh = @() val;
%         config.gh = @(p) zeros(objLen, 1);
%         config.str = strcat('%s',config.type,str);
%         return;
%     end
    
    % option 2: it is single variable
%     rhv = cellfun(@(c) strcmp(str, c), varNames);
%     if sum(rhv)==1
%         config.rhv = find(rhv);
%         config.fh = @(p) p;
%         config.gh = @(jac, ~) jac;
%         config.str = strcat('%s',config.type,'%s');
%         return;
%     end
    
    % parse the right-hand expression
    varInput = strjoin(cellfun(@(c) strcat('(',c,')?'), varNames, 'UniformOutput', false),'');
    numInput = strcat('(pi)|(\d*(\.)?\d*(e-)?(e+)?(e)?\d*)');
    format = strcat(numInput, varInput);
    [match, nomatch] = regexp(str, format, 'match', 'split');
    matchLen = length(match);
    isBasicOperator = ismember(nomatch(2:matchLen),{'-', '+', '*', '/'});
    assert(all(isBasicOperator)&&isempty(nomatch{matchLen+1}),...
        'MATLAB:linker:parseConstraint', ...
        'Unknown format.');
    
    % define rhv
    isVariable = cellfun(@(c) whichVar(c, varNames), match);
    rhv = isVariable(isVariable>0);
    assert(length(unique(rhv))==length(rhv),...
        'MATLAB:linker:parseConstraint',...
        'Each variable must be used once in each constraint.');
    config.rhv = rhv;
    
    % define constraint str
    tmp = match; tmp(isVariable>0) = {'%s'};
    config.str = strcat('%s', config.type, nomatch{1}, ...
        strjoin(tmp, nomatch(2:matchLen)));
    
    % define function handle
    tmp = match;
    for n=1:matchLen
        if isVariable(n)>0
            tmp{n} = sprintf('P(%d)', find(rhv==isVariable(n)));
        end
    end
    config.fh = str2func(strcat('@(P)',nomatch{1},...
        strjoin(tmp, nomatch(2:matchLen))));
    
    % define gradient handle
    tmpMatch = match(isVariable>0);
    tmpNomatch = nomatch(isVariable>0);
    tmpStr = '';
    for n=1:length(tmpMatch)
        tmpStr = strcat(tmpStr, tmpNomatch{n}, sprintf('jac(%d,:)', n));
    end
    if isempty(tmpStr)
        config.gh = @(~,~) [];
    else
        config.gh = str2func(strcat('@(jac,~)',tmpStr));
    end
end

function i = whichVar(thisVar, varNames)
    i = find(strcmp(thisVar, varNames));
    if isempty(i)
        i = 0;
    end
end