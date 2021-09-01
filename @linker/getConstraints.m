 function constraintList = getConstraints(obj)
    % getConstraints is a public method for the class linker
    %
    %   getConstraints(obj) returns a list of constraints applied to
    %   optimisation variables.
    %       Input Arguments: 
    %                   obj: A column array of class linker
    %      Output Arguments:
    %        constraintList: A cell-string
    
    constraintList = {};
    [~,index] = obj.getCompOrder();
    for n = index'
        switch obj(n).type
            case 'independent'
                if obj(n).bounded
                    lstr = sprintf('%s>=%1.4e', obj(n).name, obj(n).minimum);
                    ustr = sprintf('%s<=%1.4e', obj(n).name, obj(n).maximum);
                    constraintList = [constraintList; lstr; ustr];
                end
                %
            case 'dependent'
                idx = obj(n).dependence;
                varNames = {obj([n, idx]).name};
                if contains(obj(n).constraint, '>=')
                    lstr = sprintf(obj(n).constraint, varNames{:});
                    ustr = sprintf('%s<=%1.4e', obj(n).name, obj(n).maximum);
                else
                    lstr = sprintf('%s>=%1.4e', obj(n).name, obj(n).minimum);
                    ustr = sprintf(obj(n).constraint, varNames{:});
                end
                constraintList = [constraintList; lstr; ustr];
                %
            case 'constant'
                str = sprintf(obj(n).constraint, obj(n).name);
                constraintList = [constraintList; str];
                %
            case 'dummy'
                idx = obj(n).dependence;
                varNames = {obj([n, idx]).name};
                str = sprintf(obj(n).constraint, varNames{:});
                constraintList = [constraintList; str];
                %
            otherwise
                error('MATLAB:linker:getConstraints',...
                    'Unknow linker type');
        end%swith-case
    end%for
end