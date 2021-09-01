function addConstraint(obj, str)
    % addConstraint is a method for class LINKER. 
    % 
    % addConstraint(obj, str) enforce a new constraint over the object
    % variables.
    %       Input Arguments:
    %                   obj: A column linker object
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
    
    % parse input constraint
    config = parseConstraint(obj, str);
    
    % equality constraint, constant value
    if strcmp(config.type, '=') && isempty(config.rhv)
        index = config.lhv;
        C = config.fh();
        if obj(index).bounded
            assert(C<=obj(index).maximum && C>=obj(index).minimum,...
                'MATLAB:linker:addConstraint',...
                '%s must be between %1.4e and %1.4e.',...
                obj(index).name, obj(index).minimum, obj(index).maximum);
        end
        
        obj(index).type = 'constant';
        obj(index).dependence = [];
        obj(index).constraint = config.str;
        obj(index).fval = @(~,~,~,~) C;
        obj(index).fvalinv = @(~,~,~) 0;
        obj(index).fgrad = @(~,~,~) 0;
        obj(index).ubound = @(~,~,~) [];
        obj(index).lbound = @(~,~,~) [];
        obj(index).uboundgrad = @(~,~,~) [];
        obj(index).lboundgrad = @(~,~,~) [];
        obj.updateCompOrder();
        return;
    end
    
    % equality constraint, dummy variable
    if strcmp(config.type, '=') && not(isempty(config.rhv))
        index = config.lhv;
        obj(index).type = 'dummy';
        obj(index).constraint = config.str;
        obj(index).dependence = config.rhv;
        obj(index).fval = @(~,P,~,~) config.fh(P);
        obj(index).fvalinv = @(~,~,~) 0;
        obj(index).fgrad = @(jac, P) config.gh(jac, P);
        
        obj(index).ubound = @(~,~,~) [];
        obj(index).lbound = @(~,~,~) [];
        obj(index).uboundgrad = @(~,~,~) [];
        obj(index).lboundgrad = @(~,~,~) [];
        obj.updateCompOrder();
        return;
    end
    
    % set the lowerbound
    if strcmp(config.type, '>=') && isempty(config.rhv)
        index = config.lhv;
        C = config.fh();

        assert(obj(index).bounded, ...
            'MATLAB:linker:addConstraint',...
            '%s must be bounded to set a lower bound.',...
            obj(index).name);
        assert(obj(index).maximum>C, ...
            'MATLAB:linker:addConstraint',...
            'Input lower bound must be less than %1.4e.',...
            obj(index).maximum);
        
        obj(index).minimum = C;
        return;
    end
    
    % dependent parameter where lowerband is a function of other variables
    if strcmp(config.type, '>=') && not(isempty(config.rhv))
        index = config.lhv;
        assert(obj(index).bounded && strcmp(obj(index).type,'independent'), ...
            'MATLAB:linker:addConstraint',...
            '%s must be bounded and independent to enforce inequality constraint.',...
            obj(index).name);
        idx = obj.getDependent(config.rhv);
        assert(sum(index==idx)==0, 'MATLAB:linker:addConstraint',...
            'Right-hand expression must be independent of %s.',...
            obj(index).name);
        
        obj(index).type = 'dependent';
        obj(index).constraint = config.str;
        obj(index).dependence = config.rhv;
        obj(index).lbound = @(p,~,~) config.fh(p);
        obj(index).lboundgrad = @(jac, p, ~) config.gh(jac, p);
        obj.updateCompOrder();
        return;
    end
    
    % set the upperbound
    if strcmp(config.type, '<=') && isempty(config.rhv)
        index = config.lhv;
        C = config.fh();

        assert(obj(index).bounded, ...
            'MATLAB:linker:addConstraint',...
            '%s must be bounded to set a lower bound.',...
            obj(index).name);
        assert(obj(index).minimum<C, ...
            'MATLAB:linker:addConstraint',...
            'Input upper bound must be greater than %1.4e.',...
            obj(index).minimum);
        
        obj(index).maximum = C;
        return;
    end
    
    % dependent parameter where upperband is a function of other variables
    if strcmp(config.type, '<=') && not(isempty(config.rhv))
        index = config.lhv;
        assert(obj(index).bounded && strcmp(obj(index).type,'independent'), ...
            'MATLAB:linker:addConstraint',...
            '%s must be bounded and independent to enforce inequality constraint.',...
            obj(index).name);
        idx = obj.getDependent(config.rhv);
        assert(sum(index==idx)==0, 'MATLAB:linker:addConstraint',...
            'Right-hand expression must be independent of %s.',...
            obj(index).name);
        
        obj(index).type = 'dependent';
        obj(index).constraint = config.str;
        obj(index).dependence = config.rhv;
        obj(index).ubound = @(p,~,~) config.fh(p);
        obj(index).uboundgrad = @(jac, p, ~) config.gh(jac, p);
        obj.updateCompOrder();
        return;
    end
    
    %--- pi=c -------------------------------------------------------------
    % check if str is a constant constraint
    [v, ~, msg] = sscanf(str, 'p%d=%f');
    if isempty(msg) && length(v) == 2
        % input string is matched with the given format
        fixedLinkNo = v(1);
        fixedLinkVal = v(2);

        % check if fixedLink is a valid link number
        assert(fixedLinkNo>0 && fixedLinkNo<length(obj)+1, ...
               sprintf(['Input link should be an interger ', ...
                        'between 1 and %d.'], length(obj)));

        % assert that fixedLink is not an dummy
        if obj(fixedLinkNo).dummy
            error('Link %d is an dummy and cannot be fixed.', ...
                  fixedLinkNo);
        end

        % assert that fixedLink is not used in an inequality
        if obj(fixedLinkNo).linked>0
            error(['Link %d is part of an inequality ', ...
                   'constraint and cannot be fixed.'], ...
                   fixedLinkNo);
        end


        % update values                 
        obj(fixedLinkNo).type = 'constant';
        obj(fixedLinkNo).storedParamVal = fixedLinkVal;
        obj(fixedLinkNo).constant = true;
        obj(fixedLinkNo).lowerBound = -inf;
        obj(fixedLinkNo).upperBound = inf;
        obj.sortCompOrder();
        return;
    end
    
    %--- pi>=c ------------------------------------------------------------
    [v, ~, msg] = sscanf(str, 'p%d>=%f');
    if isempty(msg) && length(v) == 2
        linkNo = v(1);
        numValue = v(2);
        obj(linkNo).set('lowerBound', numValue);
        return;
    end
    
    %--- pi<=c ------------------------------------------------------------
    [v, ~, msg] = sscanf(str, 'p%d<=%f');
    if isempty(msg) && length(v) == 2
        linkNo = v(1);
        numValue = v(2);
        obj(linkNo).set('upperBound', numValue);
        return;
    end
    
    %--- pi=pj ------------------------------------------------------------
    % check if str is an equality constraint
    [v, ~, msg] = sscanf(str, 'p%d=p%d');
    if isempty(msg) && length(v) == 2
        % input string is matched with the given format
        dummyLink = v(1);
        if obj(dummyLink).linked>0 || obj(dummyLink).constant
            error(['Parameter p%d is not free and cannot ', ...
                   'be set as a dummy link.'], dummyLink);
        end

        refLink = v(2);
        if obj(refLink).dummy
            error(['Parameer p%d is dummy and cannot serve ', ...
                   'as a ref link in an equality constraint.'], refLink);
        end

        obj(dummyLink).type = 'dummy';
        obj(dummyLink).dummy    = true;
        obj(dummyLink).linked   = 0;
        obj(dummyLink).constant = false;

        obj(dummyLink).crossLinkNo = refLink;
        obj(dummyLink).fval = @(thisP, ii) thisP(ii(1));
        obj(dummyLink).gval = @(thisJac, thisP, ii) thisJac(ii(1),:);
        obj(dummyLink).upperBound = inf;
        obj(dummyLink).lowerBound = -inf;
        obj(dummyLink).initRange  = [0,1];
        obj(dummyLink).symbolicFun = @(str) strcat(str{1},'=',str{2});
        obj.sortCompOrder();
        return;
    end
    
    %--- pi=(pj+pk)*c -----------------------------------------------------
    [v, ~, msg] = sscanf(str, 'p%d=(p%d+p%d)*%f');
    if isempty(msg) && length(v) == 4
        % input string is matched with the given format
        dummyLink = v(1);
        if obj(dummyLink).linked>0 || obj(dummyLink).constant
            error(['Parameter p%d is not free and cannot ', ...
                   'be set as a dummy link.'], dummyLink);
        end

        refLink1 = v(2);
        if obj(refLink1).dummy
            error(['Link %d is an alias and cannot serve ', ...
                   'as a ref link in an equality constraint.'], refLink1);
        end

        refLink2 = v(3);
        if obj(refLink2).dummy
            error(['Link %d is dummy  and cannot serve ', ...
                   'as a ref link in analias equality constraint.'], refLink2);
        end

        c = v(4); 

        obj(dummyLink).type = 'dummy';
        obj(dummyLink).dummy = true;
        obj(dummyLink).linked   = 0;
        obj(dummyLink).constant = false;

        obj(dummyLink).crossLinkNo = [refLink1, refLink2];
        obj(dummyLink).fval = @(thisP, ii) (thisP(ii(1))+thisP(ii(2)))*c;
        obj(dummyLink).gval = @(thisJac, thisP, ii) (thisJac(ii(1),:)+thisJac(ii(2)))*c;
        obj(dummyLink).upperBound = inf;
        obj(dummyLink).lowerBound = -inf;
        obj(dummyLink).initRange  = [0,1];
        obj(dummyLink).symbolicFun = @(str) strcat(str{1},'=(',str{2}, '+', str{3}, ')*', num2str(c));
        obj.sortCompOrder();
        return;
    end
    
    %--- pi= pj/(pk+pl) ---------------------------------------------------
    [v, ~, msg] = sscanf(str, 'p%d=p%d/(p%d+p%d)');
    if isempty(msg) && length(v) == 4
        % input string is matched with the given format
        dummyLink = v(1);
        if obj(dummyLink).linked>0 || obj(dummyLink).constant
            error(['Parameter p%d is not free and cannot ', ...
                   'be set as a dummy link.'], dummyLink);
        end

        refLink1 = v(2);
        if obj(refLink1).dummy
            error(['Link %d is an alias and cannot serve ', ...
                   'as a ref link in an equality constraint.'], refLink1);
        end

        refLink2 = v(3);
        if obj(refLink2).dummy
            error(['Link %d is dummy  and cannot serve ', ...
                   'as a ref link in analias equality constraint.'], refLink2);
        end

        refLink3 = v(4);
        if obj(refLink3).dummy
            error(['Link %d is dummy  and cannot serve ', ...
                   'as a ref link in analias equality constraint.'], refLink3);
        end


        obj(dummyLink).type = 'dummy';
        obj(dummyLink).dummy = true;
        obj(dummyLink).linked   = 0;
        obj(dummyLink).constant = false;

        obj(dummyLink).crossLinkNo = [refLink1, refLink2, refLink3];
        obj(dummyLink).fval = @(thisP, ii) thisP(ii(1))/(thisP(ii(2))+thisP(ii(3)));
        obj(dummyLink).gval = @(thisJac, thisP, ii) ...
            (thisJac(ii(1),:)*(thisP(ii(2))+thisP(ii(3)))-(thisJac(ii(2),:)+thisJac(ii(3)))*thisP(ii(1)))/(thisP(ii(2))+thisP(ii(3)))^2;
        obj(dummyLink).upperBound = inf;
        obj(dummyLink).lowerBound = -inf;
        obj(dummyLink).initRange  = [0,1];
        obj.sortCompOrder();
        return;
    end
    
    %--- pi=(1-pj)*c -----------------------------------------------------
    [v, ~, msg] = sscanf(str, 'p%d=(1-p%d)*%f');
    if isempty(msg) && length(v) == 3
        % input string is matched with the given format
        dummyLink = v(1);
        if obj(dummyLink).linked>0 || obj(dummyLink).constant
            error(['Parameter p%d is not free and cannot ', ...
                   'be set as a dummy link.'], dummyLink);
        end

        refLink1 = v(2);
        if obj(refLink1).dummy
            error(['Link %d is an alias and cannot serve ', ...
                   'as a ref link in an equality constraint.'], refLink1);
        end

        c = v(3); 

        obj(dummyLink).type = 'dummy';
        obj(dummyLink).dummy = true;
        obj(dummyLink).linked   = 0;
        obj(dummyLink).constant = false;

        obj(dummyLink).crossLinkNo = [refLink1];
        obj(dummyLink).fval = @(thisP, ii) (1-thisP(ii(1)))*c;
        obj(dummyLink).gval = @(thisJac, thisP, ii) -thisJac(ii(1),:)*c;
        obj(dummyLink).upperBound = inf;
        obj(dummyLink).lowerBound = -inf;
        obj(dummyLink).initRange  = [0,1];
        obj.sortCompOrder();
        return;
    end
    
    %--- pi=pj*c ----------------------------------------------------------
    [v, ~, msg] = sscanf(str, 'p%d=p%d*%f');
    if isempty(msg) && length(v) == 3
        % input string is matched with the given format
        dummyLink = v(1);
        if obj(dummyLink).linked>0 || obj(dummyLink).constant
            error(['Parameter p%d is not free and cannot ', ...
                   'be set as a dummy link.'], dummyLink);
        end

        refLink1 = v(2);
        if obj(refLink1).dummy
            error(['Link %d is an alias and cannot serve ', ...
                   'as a ref link in an equality constraint.'], refLink1);
        end

        c = v(3); 

        obj(dummyLink).type = 'dummy';
        obj(dummyLink).dummy = true;
        obj(dummyLink).linked   = 0;
        obj(dummyLink).constant = false;

        obj(dummyLink).crossLinkNo = [refLink1];
        obj(dummyLink).fval = @(thisP, ii) thisP(ii(1))*c;
        obj(dummyLink).gval = @(thisJac, thisP, ii) thisJac(ii(1),:)*c;
        obj(dummyLink).upperBound = inf;
        obj(dummyLink).lowerBound = -inf;
        obj(dummyLink).initRange  = [0,1];
        obj.sortCompOrder();
        return;
    end
    
    %--- pi=(1-pj-pk)*c ---------------------------------------------------
    [v, ~, msg] = sscanf(str, 'p%d=(1-p%d-p%d)*%f');
    if isempty(msg) && length(v) == 4
        % input string is matched with the given format
        dummyLink = v(1);
        if obj(dummyLink).linked>0 || obj(dummyLink).constant
            error(['Parameter p%d is not free and cannot ', ...
                   'be set as a dummy link.'], dummyLink);
        end

        refLink1 = v(2);
        if obj(refLink1).dummy
            error(['Link %d is an alias and cannot serve ', ...
                   'as a ref link in an equality constraint.'], refLink1);
        end

        refLink2 = v(3);
        if obj(refLink2).dummy
            error(['Link %d is an alias and cannot serve ', ...
                   'as a ref link in an equality constraint.'], refLink2);
        end

        c = v(4); 

        obj(dummyLink).type = 'dummy';
        obj(dummyLink).dummy = true;
        obj(dummyLink).linked   = 0;
        obj(dummyLink).constant = false;

        obj(dummyLink).crossLinkNo = [refLink1, refLink2];
        obj(dummyLink).fval = @(thisP, ii) (1-thisP(ii(1))-thisP(ii(2)))*c;
        obj(dummyLink).gval = @(thisJac, thisP, ii) -(thisJac(ii(1),:)+thisJac(ii(2),:))*c;
        obj(dummyLink).upperBound = inf;
        obj(dummyLink).lowerBound = -inf;
        obj(dummyLink).initRange  = [0,1];
        obj.sortCompOrder();
        return;
    end
    
    %--- pi=pj*pk ---------------------------------------------------------
    [v, ~, msg] = sscanf(str, 'p%d=p%d*p%d');
    if isempty(msg) && length(v) == 3
        % input string is matched with the given format
        dummyLink = v(1);
        if obj(dummyLink).linked>0 || obj(dummyLink).constant
            error(['Parameter p%d is not free and cannot ', ...
                   'be set as a dummy link.'], dummyLink);
        end

        refLink1 = v(2);
        if obj(refLink1).dummy
            error(['Link %d is an alias and cannot serve ', ...
                   'as a ref link in an equality constraint.'], refLink1);
        end

        refLink2 = v(3);
        if obj(refLink2).dummy
            error(['Link %d is dummy  and cannot serve ', ...
                   'as a ref link in analias equality constraint.'], refLink2);
        end

        obj(dummyLink).type = 'dummy';
        obj(dummyLink).dummy = true;
        obj(dummyLink).linked   = 0;
        obj(dummyLink).constant = false;

        obj(dummyLink).crossLinkNo = [refLink1, refLink2];
        obj(dummyLink).fval = @(thisP, ii) thisP(ii(1))*thisP(ii(2));
        obj(dummyLink).gval = @(thisJac, thisP, ii) thisP(ii(1))*thisJac(ii(2),:) + thisP(ii(2))*thisJac(ii(1));
        obj(dummyLink).upperBound = inf;
        obj(dummyLink).lowerBound = -inf;
        obj(dummyLink).initRange  = [0,1];
        obj.sortCompOrder();
        return;
    end
    
    %--- pi=pj*(1-pk) -----------------------------------------------------
    [v, ~, msg] = sscanf(str, 'p%d=p%d*(1-p%d)');
    if isempty(msg) && length(v) == 3
        % input string is matched with the given format
        dummyLink = v(1);
        if obj(dummyLink).linked>0 || obj(dummyLink).constant
            error(['Parameter p%d is not free and cannot ', ...
                   'be set as a dummy link.'], dummyLink);
        end

        refLink1 = v(2);
        if obj(refLink1).dummy
            error(['Link %d is an alias and cannot serve ', ...
                   'as a ref link in an equality constraint.'], refLink1);
        end

        refLink2 = v(3);
        if obj(refLink2).dummy
            error(['Link %d is dummy  and cannot serve ', ...
                   'as a ref link in analias equality constraint.'], refLink2);
        end

        obj(dummyLink).type = 'dummy';
        obj(dummyLink).dummy = true;
        obj(dummyLink).linked   = 0;
        obj(dummyLink).constant = false;

        obj(dummyLink).crossLinkNo = [refLink1, refLink2];
        obj(dummyLink).fval = @(thisP, ii) thisP(ii(1))*(1-thisP(ii(2)));
        obj(dummyLink).gval = @(thisJac, thisP, ii) -thisP(ii(1))*thisJac(ii(2),:) + (1-thisP(ii(2)))*thisJac(ii(1));
        obj(dummyLink).upperBound = inf;
        obj(dummyLink).lowerBound = -inf;
        obj(dummyLink).initRange  = [0,1];
        obj(dummyLink).symbolicFun = @(str) strcat(str{1},'=',str{2},'*(1-',str{3},')');
        obj.sortCompOrder();
        return;
    end
    
    %--- pi=pj*(1-pk-pl) --------------------------------------------------
    [v, ~, msg] = sscanf(str, 'p%d=p%d*(1-p%d-p%d)');
    if isempty(msg) && length(v) == 4
        % input string is matched with the given format
        dummyLink = v(1);
        if obj(dummyLink).linked>0 || obj(dummyLink).constant
            error(['Parameter p%d is not free and cannot ', ...
                   'be set as a dummy link.'], dummyLink);
        end

        refLink1 = v(2);
        if obj(refLink1).dummy
            error(['Link %d is an alias and cannot serve ', ...
                   'as a ref link in an equality constraint.'], refLink1);
        end

        refLink2 = v(3);
        if obj(refLink2).dummy
            error(['Link %d is dummy  and cannot serve ', ...
                   'as a ref link in analias equality constraint.'], refLink2);
        end

        refLink3 = v(4);
        if obj(refLink3).dummy
            error(['Link %d is dummy  and cannot serve ', ...
                   'as a ref link in analias equality constraint.'], refLink3);
        end

        obj(dummyLink).type = 'dummy';
        obj(dummyLink).dummy = true;
        obj(dummyLink).linked   = 0;
        obj(dummyLink).constant = false;

        obj(dummyLink).crossLinkNo = [refLink1, refLink2, refLink3];
        obj(dummyLink).fval = @(thisP, ii) thisP(ii(1))*(1-thisP(ii(2))-thisP(ii(3)));
        obj(dummyLink).gval = @(thisJac, thisP, ii) -thisP(ii(1))*(thisJac(ii(2),:)+thisJac(ii(3),:)) + (1-thisP(ii(2))-thisP(ii(3)))*thisJac(ii(1));
        obj(dummyLink).upperBound = inf;
        obj(dummyLink).lowerBound = -inf;
        obj(dummyLink).initRange  = [0,1];
        obj.sortCompOrder();
        return;
    end
    
    %--- pi>=pj -----------------------------------------------------------
    % check if str is an inequality constraint
    [v, ~, msg] = sscanf(str, 'p%d>=p%d');
    if isempty(msg) && length(v) == 2
        % input string is matched with the given format
        greaterLink = v(1);
        smallerLink = v(2);

        smallerTypeIsNotCos = ~strcmp(obj(smallerLink).type, 'cos');
        if smallerTypeIsNotCos
            error(['For ineuality constraints ">=", left ', ...
                   'argument should be of tyep "cos".']);
        end

        if obj(greaterLink).dummy || obj(smallerLink).dummy
            error('Input links cannot be dummy.');
        end

%                 if strcmp(obj(smallerLink).type, 'boundedCos')
%                     error(['Link %d cannot serve as a smaller link ', ...
%                            'in an inequality constraint']);
%                 end

        obj(smallerLink).crossLinkNo = greaterLink;
        obj(smallerLink).linked = 1;
        obj(smallerLink).dummy = false;
        obj(smallerLink).constant = false;

        if obj(greaterLink).compOrder > obj(smallerLink).compOrder
            tmp = obj(greaterLink).compOrder;
            obj(greaterLink).compOrder = obj(smallerLink).compOrder;
            obj(smallerLink).compOrder = tmp;
        end
        return;
    end
    
    %--- pi+pj=c ----------------------------------------------------------
     % check if str is a cumsum2 constraint
    [v, ~, msg] = sscanf(str, 'p%d+p%d=%f');
    if isempty(msg) && length(v) == 3
        % input string is matched with the given format
        link1No = v(1);
        link2No = v(2);
        C = v(3);

        obj.isValidInputForCumsum2(link1No, link2No, C);

        obj(link1No).upperBound  = C - obj(link2No).lowerBound;
        obj(link2No).upperBound  = C - obj(link1No).lowerBound;


        obj(link2No).type = 'dummy';
        obj(link2No).dummy = true;
        obj(link2No).constant = false;
        obj(link2No).linked = 0;
        obj(link2No).crossLinkNo = link1No;
        obj(link2No).fval = @(thisP, ii) C-thisP(ii(1));
        obj(link2No).gval = @(thisJac, thisP, ii) -thisJac(ii(1),:); 
        obj(link2No).symbolicFun = @(str) strcat(str{1},'+',str{2},'=',num2str(C));
        return;
    end
    
    %--- pi+pj+pk=c -------------------------------------------------------
    % check if str is a cumsum3 constraint
    [v, ~, msg] = sscanf(str, 'p%d+p%d+p%d=%f');
    if isempty(msg) && length(v) == 4
        % input string is matched with the given format
        link1No = v(1);
        link2No = v(2);
        link3No = v(3);
        C = v(4);

        obj.isValidInputForCumsum3(link1No, link2No, link3No, C);

        obj(link1No).upperBound  = C - obj(link2No).lowerBound - obj(link3No).lowerBound;
        obj(link2No).upperBound  = C - obj(link3No).lowerBound;
        %obj(link3No).upperBound  = C - obj(link1No).lowerBound - obj(link2No).lowerBound;


        obj(link2No).linked = 2;
        obj(link2No).crossLinkNo = link1No;

        if obj(link1No).compOrder > obj(link2No).compOrder
            tmp = obj(link1No).compOrder;
            obj(link1No).compOrder = obj(link2No).compOrder;
            obj(link2No).compOrder = tmp;
        end

        obj(link3No).type = 'dummy';
        obj(link3No).dummy = true;
        obj(link3No).crossLinkNo = [link1No, link2No];
        obj(link3No).fval = @(thisP, ii) C-thisP(ii(1))-thisP(ii(2));
        obj(link3No).gval = @(thisJac, thisP, ii) -thisJac(ii(1),:)-thisJac(ii(2),:);
        obj(link3No).symbolicFun = @(str) strcat(str{1},'+',str{2},'+',str{3},'=',num2str(C));
        obj.sortCompOrder();
        return;
    end

    % the constraint has not been recognised yet!
    error('Constraint "%s" is not recognised!', str);
end % of addConstraint

function [flag, thisLinkNo, crossLinkNo, fh] = isEqual(str)
    % isEqual check if the input constraint is of type equality
    % pi = fh(pj, pk, ...)
    
    % initialise with empty vector
    thisLinkNo = [];
    crossLinkNo = [];
    fh = [];
    
    % split the left-hand and the right-hand sides of equation
    [v, thisOperator] = strsplit(str, {'=', '>', '<'});
    
    % identify if the constraint is of equality type; otherwise return
    % false
    isOperatorOk = length(thisOperator)==1 && strcmp(thisOperator{1}, '=');
    if ~isOperatorOk
        flag = false;
        return;
    else
        leftSide = v{1};
        rightSide = v{2};
    end
    
    % check if the link number on the left side is valid    
    [thisLinkNo, isLeftSideOk] = str2linkNo(leftSide);
    if ~isLeftSideOk
        flag = false;
        return;
    end
    
    % Extract crossLinkNo from the right side of equation
    [variables, operators] = strsplit(rightSide, {')','(','+','-','*','/'});
    nVariables = length(variables);
    %
    for idx = 1:nVariables
        counter = 0;
        % check if it is a constant number then skip
        if isnan(str2double(variables{idx}))
            [refLink, flag] = str2linkNo(variables{idx});
            crossLinkNo = [crossLinkNo, refLink];
            % if conversion is not successful raise error
            if not(flag)
                error('MATlAB:linker:addConstraint',...
                      "Constraint '%s' is not recognised.", str);
            end
            % replace variable name with x(idx) for use in str2func
            counter = counter + 1;
            variables{idx} = sprintf('x(%d)', counter);
        end
    end
    
    % create function handles
    fh = str2func(strcat('@(x)', strjoin(variables, operators)));

end

function [linkNo, flag] = str2linkNo(str)
    [linkNo, ~, msg] = sscanf(str, 'p%d');
    flag = isempty(msg) && isscalar(linkNo);
end

