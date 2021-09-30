function idx = select(obj, str)
    % select is public method for class scheme
    
    assert(ischar(str), 'MATLAB:scheme:invalidInputArgument',...
        'Input must be of type char.');
    
    strTree = generateTree(str);
    
    % convert strings into boolean variables in leaves
    id = strTree.findleaves;
    for i = id
        this_str = strTree.get(i);
        idx = evaluateSelectionFilter(obj, this_str);
        strTree = strTree.set(i, idx);
    end
    
        % merge compartments together from bottom to top
    while strTree.nnodes > 1
        leaves = arrayfun(@(i) strTree.isleaf(i), 1:strTree.nnodes);
        
        for i = find(leaves)
            parentId = strTree.getparent(i);
            siblings = strTree.getsiblings(i);
            parentLable = strTree.get(parentId);
            if strcmp(parentLable, 'C')
                assert(siblings==i, 'MATLAB:scheme:unknownError',...
                    'Unknown error with the strTree structure.');
                idx = strTree.get(i);
                strTree = strTree.set(parentId, idx);
                strTree = strTree.remove(i);
                break;
            elseif strcmp(parentLable, 'AND')
                assert(length(siblings)==2, 'MATLAB:scheme:unknownError',...
                    'Unknown error with the strTree structure.');
                idx = strTree.Node(siblings);
                if islogical(idx{1}) && islogical(idx{2})
                    strTree = strTree.set(parentId, idx{1}&idx{2});
                    % remove nodes
                    siblingTags = strTree.tag(siblings);
                    for j = 1:2
                        thisNode = find(cellfun(@(c) strcmp(c, siblingTags{j}), strTree.tag));
                        strTree = strTree.removenode(thisNode);
                    end
                    break;
                end
            elseif strcmp(parentLable, 'OR')
                assert(length(siblings)==2, 'MATLAB:scheme:unknownError',...
                    'Unknown error with the strTree structure.');
                idx = strTree.Node(siblings);
                if islogical(idx{1}) && islogical(idx{2})
                    strTree = strTree.set(parentId, idx{1}|idx{2});
                    % remove nodes
                    siblingTags = strTree.tag(siblings);
                    for j = 1:2
                        thisNode = find(cellfun(@(c) strcmp(c, siblingTags{j}), strTree.tag));
                        strTree = strTree.removenode(thisNode);
                    end
                    break;
                end
            end
        end
    end
       
    idx = strTree.get(1);
end

function idx = evaluateSelectionFilter(obj, str)
    % return boolean variable
    
    [match, nomatch] = regexp(str, '(>=)?|(<=)?|(>)?|(<)?|(=)?', 'match', 'split');
    
    assert(isscalar(match) && all(size(nomatch)==[1, 2]),...
        'MATLAB:scheme:invalidInputArgument',...
        'Input "%s" is not valid as a selection criteria.', str);
    
    lhv = obj.get(nomatch{1});
    rhv = str2num(nomatch{2});
    operator = match{1};
    
    switch operator
        case '='
            assert(isvector(rhv)&&~isempty(rhv)&&isnumeric(rhv),...
                'MATLAB:scheme:invalidInputArgument',...
                'Input %s is not valid.', str);
            idx = ismember(lhv, rhv);
        case '>='
            assert(isscalar(rhv)&&~isempty(rhv)&&isnumeric(rhv),...
                'MATLAB:scheme:invalidInputArgument',...
                'Input %s is not valid.', str);
            idx = lhv>=rhv;
        case '<='
            assert(isscalar(rhv)&&~isempty(rhv)&&isnumeric(rhv),...
                'MATLAB:scheme:invalidInputArgument',...
                'Input %s is not valid.', str);
            idx = lhv<=rhv;    
        case '>'
            assert(isscalar(rhv)&&~isempty(rhv)&&isnumeric(rhv),...
                'MATLAB:scheme:invalidInputArgument',...
                'Input %s is not valid.', str);
            idx = lhv>rhv; 
        case '<'
            assert(isscalar(rhv)&&~isempty(rhv)&&isnumeric(rhv),...
                'MATLAB:scheme:invalidInputArgument',...
                'Input %s is not valid.', str);
            idx = lhv<rhv;
        otherwise
            error('MATLAB:scheme:invalidInputArgument',...
                'Input operator %s is not valid in %s.', operator, str);
    end
end

function strTree = generateTree(str)
    % Generate a tree structure to represent the required hierarchy to 
    % compute logical operations AND and OR in the input string
    
    % string pre-processing
    str(isspace(str)) = []; % remove white space
    str = strrep(str, 'AND', '&');
    str = strrep(str, 'OR', '|');
    
    % first node
    strTree = tree('C');
    parentId = 1;
    
    strLen = length(str);
    for i = 1:strLen
        switch str(i)
            case '('
                if strcmp(strTree.get(parentId), 'C')
                    strTree = strTree.addnode(parentId, 'C');
                    parentId = strTree.nnodes;
                else
                    parentId = strTree.nnodes;
                    strTree = strTree.set(parentId, 'C');
                end
                
            case ')'
                parentId = strTree.Parent(parentId);
                
            case '&'
                nodeLable = strTree.get(parentId);
                if strcmp(nodeLable, 'C')
                    strTree = strTree.set(parentId, 'AND');
                    % add a new node
                    strTree = strTree.addnode(parentId, '');
                elseif ismember(nodeLable, {'AND', 'OR'})
                    % insert a new node
                    child = parentId;
                    parent = strTree.Parent(child);
                    [strTree, parentId] = strTree.insertNode(parent, child, 'AND');
                    % add a new node
                    strTree = strTree.addnode(parentId, '');
                else
                    error('MATLAB:scheme:unexpected',...
                        'This is an unknown error.');
                end
                
            case '|'
                nodeLable = strTree.get(parentId);
                if strcmp(nodeLable, 'C')
                    strTree = strTree.set(parentId,'OR');
                    % add a new node
                    strTree = strTree.addnode(parentId, '');
                elseif ismember(nodeLable, {'AND', 'OR'})
                    % insert a new node
                    child = parentId;
                    parent = strTree.Parent(child);
                    [strTree, parentId] = strTree.insertNode(parent, child, 'OR');
                    % add a new node
                    strTree = strTree.addnode(parentId, '');
                else
                    error('MATLAB:scheme:unexpected',...
                        'This is an unknown error.');
                end
                
            otherwise
                % it is an experession
                lastNode = strTree.nnodes;
                nodeLable = strTree.get(lastNode);
                if ismember(nodeLable, {'C', 'AND', 'OR'})
                    % add a new node
                    strTree = strTree.addnode(parentId, str(i));
                else
                    % update the node lable
                    nodeLable = strcat(nodeLable, str(i));
                    strTree = strTree.set(lastNode, nodeLable);
                end
        end% switch-case
    end% for loop
end% generateTree(str)
