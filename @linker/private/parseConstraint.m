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
    
    %% parse the right-hand expression
    % include a dummy zero if str starts with + or -
    if strcmp(str(1), '+') || strcmp(str(1), '-')
        str = append('0', str);
    end
    % split the construction into additive or substractive expressions
    [operators, expressions] = regexp(str, '(?<![\d*]e)+|(?<![\d*]e)-',...
                                      'match', 'split');
    
    rhv = [];
    str_expr = '';
    str_fh = '';
    str_gh = '';
    nOperators = length(operators);
    for i=1:length(expressions)
        [rhv, str_expr, str_fh, str_gh] = ...
            readExpression(expressions{i}, varNames, rhv, str_expr, str_fh, str_gh);
        if i<=nOperators
            str_expr = append(str_expr,operators{i});
            str_fh   = append(str_fh  ,operators{i});
            str_gh   = append(str_gh  ,operators{i});
        end
    end
    
    % define dependent variables
    config.rhv = rhv;
    % define constraint str
    config.str = append('%s', config.type, str_expr);
    % function handle 
    config.fh = str2func(append('@(P)',str_fh));
    % function handle to gradient
    config.gh = str2func(append('@(jac,P)',str_gh));
    
%     varInput = strjoin(cellfun(@(c) strcat('(',c,')?'), varNames, 'UniformOutput', false),'');
%     numInput = strcat('(pi)|(\d*(\.)?\d*(e-)?(e+)?(e)?\d*)');
%     format = strcat(numInput, varInput);
%     [match, nomatch] = regexp(str, format, 'match', 'split');
%     matchLen = length(match);
%     isBasicOperator = ismember(nomatch(2:matchLen),{'-', '+', '*', '/'});
%     assert(all(isBasicOperator)&&isempty(nomatch{matchLen+1}),...
%         'MATLAB:linker:parseConstraint', ...
%         'Unknown format.');
    
    % define rhv
%     isVariable = cellfun(@(c) whichVar(c, varNames), match);
%     rhv = isVariable(isVariable>0);
%     assert(length(unique(rhv))==length(rhv),...
%         'MATLAB:linker:parseConstraint',...
%         'Each variable must be used once in each constraint.');
%     config.rhv = rhv;
    
    % define constraint str
%     tmp = match; tmp(isVariable>0) = {'%s'};
%     config.str = strcat('%s', config.type, nomatch{1}, ...
%         strjoin(tmp, nomatch(2:matchLen)));
    
    % define function handle
%     tmp = match;
%     for n=1:matchLen
%         if isVariable(n)>0
%             tmp{n} = sprintf('P(%d)', find(rhv==isVariable(n)));
%         end
%     end
%     config.fh = str2func(strcat('@(P)',nomatch{1},...
%         strjoin(tmp, nomatch(2:matchLen))));
%     
%     % define gradient handle
%     tmpMatch = match(isVariable>0);
%     tmpNomatch = nomatch(isVariable>0);
%     tmpStr = '';
%     for n=1:length(tmpMatch)
%         tmpStr = strcat(tmpStr, tmpNomatch{n}, sprintf('jac(%d,:)', n));
%     end
%     if isempty(tmpStr)
%         config.gh = @(~,~) [];
%     else
%         config.gh = str2func(strcat('@(jac,~)',tmpStr));
%     end
end

function i = whichVar(thisVar, varNames)
    i = find(strcmp(thisVar, varNames));
    if isempty(i)
        i = 0;
    end
end

function [rhv, str_expr, str_fh, str_gh] = readExpression(str, varNames, rhv, str_expr, str_fh, str_gh)
    % internal function to read expressions
    % input experssion must be of types below:
    %       1) numeric input
    %       2) single parameters
    %       3) multiplication of a parameter by a numeric value
    %       4) multiplication a parameer by another one
    %       5) multiplication two paramters and a number
    
    [operators, expressions] = regexp(str, '*', 'match', 'split');
    nOperators = length(operators);
    
    switch nOperators
        case 0
            % this is single input or parameter
            varId = whichVar(expressions{1}, varNames);
            isNumber = ~isempty(str2num(expressions{1}));
            if varId>0
                rhv = [rhv, varId];
                n = length(rhv);
                str_expr = append(str_expr, '%s');
                str_fh = append(str_fh, sprintf('P(%d)', n));
                str_gh = append(str_gh, sprintf('jac(%d,:)', n));
            elseif isNumber
                str_expr = append(str_expr, expressions{1});
                str_fh = append(str_fh, expressions{1});
                str_gh = append(str_gh, '0');
            else
                error('MATLAB:linker:invalidInput',...
                'Uknown input constraint format:%s.', str);
            end
        case 1
            % two expressions: the first one could be either a number or
            % parameterthe but the second one must be a parameter 
            % 
            varId1 = whichVar(expressions{1}, varNames);
            isNumber = ~isempty(str2num(expressions{1}));
            varId2 = whichVar(expressions{2}, varNames);
            assert(varId2>0, 'MATLAB:linker:invalidInput',...
                'Uknown input constraint format:%s.', str);
            if varId1>0
                rhv = [rhv, varId1, varId2];
                n = length(rhv);
                str_expr = append(str_expr, '%s*%s');
                str_fh = append(str_fh, sprintf('P(%d)*P(%d)', n-1, n));
                str_gh = append(str_gh, sprintf('P(%d)*jac(%d,:)+P(%d)*jac(%d,:)', n-1,n,n,n-1));
            elseif isNumber
                rhv = [rhv, varId2];
                n = length(rhv);
                str_expr = append(str_expr, expressions{1},'*%s');
                str_fh = append(str_fh, expressions{1}, sprintf('*P(%d)', n));
                str_gh = append(str_gh, expressions{1}, sprintf('*jac(%d,:)', n));
            else
                error('MATLAB:linker:invalidInput',...
                'Uknown input constraint format:%s.', str);
            end
        case 2
            % three expressions: the first one must be a number and the
            % second and the third ones must be parameters
            isNumber = ~isempty(str2num(expressions{1}));
            varId1 = whichVar(expressions{2}, varNames);
            varId2 = whichVar(expressions{3}, varNames);
            assert(varId1>0 && varId2>0 && isNumber,...
                'MATLAB:linker:invalidInput',...
                'Uknown input constraint format:%s.', str);  
            rhv = [rhv, varId1, varId2];
            n = length(rhv);
            str_expr = append(str_expr, expressions{1}, '*%s*%s');
            str_fh = append(str_fh, expressions{1}, sprintf('*P(%d)*P(%d)', n-1, n));
            str_gh = append(str_gh, ...
                expressions{1}, sprintf('*P(%d)*jac(%d,:)+', n-1,n),...
                expressions{1}, sprintf('*P(%d)*jac(%d,:)', n, n-1));
        otherwise
            error('MATLAB:linker:invalidInput',...
                'Uknown input constraint format:%s.', str);
    end    
end