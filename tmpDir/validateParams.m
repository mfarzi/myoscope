function validateParams(obj, p)
% validateParams is a method for class COMPARTMENT
%
%   validateParams(obj) assert if input parameter p comply with model 
%   constraints.
%

if ~all(size(p)==[obj.nParams, 1])
    error('MATLAB:compartment:validateParams',...
          'Input must be a column vector of size %d.',...
          obj.nParams());
end

if ~isnumeric(p)
    error('MATLAB:compartment:validateParams',...
          'Input of type %s is not permitted.',...
          class(p));
end

% check consistency with constrains
x = obj.links.invLink(p);
thisP = obj.links.link(x);
err = abs(p - thisP);
idx = find(err>eps);
msg = [];
if ~isempty(idx)
    for i = idx'
        msg = strcat(msg, sprintf(['\nInput parameter %s',...
              ' does not comply with model constraints.'],...
              obj.links(i).getName{1}));
    end
    error('MATLAB:compartment:validateParams', msg);
end
end