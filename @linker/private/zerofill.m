function xp = zerofill(obj, x)
    % zerofill is a private method for class LINKER
    %
    % zerofill(obj, x) extend the input x by adding passive variables with
    % deafualt value zero.
    %
    N = size(obj, 1);
    xp = zeros(N,1);
    isActive = strcmp({obj.type}', 'dependent')|...
               strcmp({obj.type}', 'independent');
    xp(isActive) = x; 
end