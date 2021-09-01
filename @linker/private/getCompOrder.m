function [order, index] = getCompOrder(obj)
    % getCompOrder is a private method for class LINKER
    %
    % getCompOrder returns the order to compute map or invmap accounting
    % for the dependence between optimisation variables.
    %       Input argument:
    %                  obj: An array of LINKER class of size N.
    %
    %      Output Argument:
    %                order: Integer vector of size Nx1. For each array 
    %                       index obj[i], returns its queue order to
    %                       compute the map between x and p.
    %                index: Integer vector of size Nx1. Returns array 
    %                       indices in order to compute map between x and p  
    
    % read comp order
    order = arrayfun(@(thisObj) thisObj.compOrder, obj);
    % sort comp order
    [~,index] = sort(order);
end