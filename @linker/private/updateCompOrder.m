 function updateCompOrder(obj)
    % updateCompOrder is a private method for the class linker
    % 
    %   updateCompOrder(obj) sort the variables based on their dependence
    %   on other variables to allow sequential computation of variables.
    %
    
    validateattributes(obj, {'linker'}, {'column'});
    assert(all([obj.varNum]==size(obj,1)), ...
        'MATLAB:linker:updateCompOrder',...
        "'varNum' must be the same with the length of input object.");
    
    % assess dependent variables
    N = size(obj,1);
    parrents = cell(1,N);
    for n = 1:N
        parrents{n} = obj.getDependent(obj(n).dependence);
    end
    
    thisCompOrder = [];
    while(length(thisCompOrder)<N)
        for n=1:N
            parrents{n} = removeIndex(parrents{n}, thisCompOrder);
        end
        idx = cellfun(@(c) isempty(c), parrents);
        idx = removeIndex(find(idx), thisCompOrder);
        thisCompOrder = [thisCompOrder, idx];
    end
    
    % update compuation order
    for n=1:N
        obj(n).compOrder = find(thisCompOrder==n);
    end
    
    % assert all compOrders are between 1 and N
    assert(all([obj.compOrder]>=1 & [obj.compOrder]<=N), ...
        'MATLAB:linker:vertcat',...
        'Private property compOrder has bug!');
 end    

 function list1 = removeIndex(list1, list2)
    mutualMembers = ismember(list1, list2);
    list1(mutualMembers) = [];
 end