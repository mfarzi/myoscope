function model = str2model(modelName)
    if isempty(modelName)
        model = [];
        return;
    end
    
    modelTree = generateTree(modelName);
    
    % convert strings into compartment models in leaves
    id = modelTree.findleaves;
    for i = id
        thisModelName = modelTree.get(i);
        thisModel = str2compartment(thisModelName);
        modelTree = modelTree.set(i, thisModel);
    end
    
    % merge compartments together from bottom to top
    while modelTree.depth > 0
        isVoid = cellfun(@(c) strcmp(c, 'C'), modelTree.Node);
        %
        for i = 1:modelTree.nnodes
            if modelTree.isleaf(i)
                parentId = modelTree.getparent(i);
                hasVoidParent = isVoid(parentId);
                siblings = modelTree.getsiblings(i);
                hasVoidSibling = any(isVoid(siblings));
                
                if hasVoidParent && not(hasVoidSibling)
                    % merge all siblings
                    compartments = modelTree.Node(siblings);
                    thisModel = multicompartment(compartments{:});
                    modelTree = modelTree.set(parentId, thisModel);
                    
                    % remove nodes
                    siblingTags = modelTree.tag(siblings);
                    for j = 1:length(siblings)
                        thisNode = find(cellfun(@(c) strcmp(c, siblingTags{j}), modelTree.tag));
                        modelTree = modelTree.removenode(thisNode);
                    end
                    break;
                end
            end
        end % of for
    end
    
    
    model = modelTree.get(1);
    
%     if isa(modelTree.get(1), 'multicompartment')
%         model = modelTree.get(1);
%     else
%         model = multicompartment(modelTree.get(1));
%     end
    %\\
end % of str2model


function modelTree = generateTree(str)
    % Generate a tree structure to represent the required hierarchy to 
    % build the model 
    modelTree = tree('C');
    parentId = 1;
    
    strLen = length(str);
    for i = 1:strLen
        switch str(i)
            case '['
                if strcmp(str(i+1), '[')
                    s = 'C';
                else
                    s = '';
                end
                lastNode = modelTree.nnodes;
                if modelTree.get(lastNode) == 'C'
                    parentId = lastNode;
                end
                modelTree = modelTree.addnode(parentId, s); 
                
            case '-'
                if strcmp(str(i+1), '[')
                    s = 'C';
                else
                    s = '';
                end
                modelTree = modelTree.addnode(parentId, s);
                
            case ']'
                parentId = modelTree.getparent(parentId);
                
            otherwise
                lastNode = modelTree.nnodes;
                name = modelTree.get(lastNode);
                if strcmp(name, 'C')
                    name = str(i);
                else
                    name = strcat(name, str(i));
                end
                modelTree = modelTree.set(lastNode, name);
        end
    end
end

  
function model = str2compartment(modelName) 
    % transoform simple model    
    switch modelName
        case 'tensor'
            model = tensor();
        case 'zeppelin'
            model = zeppelin();
        case 'pancake'
            model = pancake();    
        case 'ball'
            model = ball();
        case 'stick'
            model = stick();    
        case 'cylinder'
            model = cylinder();
        case 'cylinderECS'
            model = cylinderECS();     
        case 'cylinderGDR'
            model = cylinderGDR();
        case 'cylinderBDA'
            model = cylinderBDA();  
        case 'zeppelinBDA'
            model = zeppelinBDA();
        case 'stickBDA'
            model = stickBDA();
        case 'zeppelinCBD'
            model = zeppelinCBD();
        otherwise
            error('model %s is not recognised', modelName);
    end            
end