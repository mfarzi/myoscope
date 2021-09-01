function n = getParamsNum(obj, mode)
% getParamsNum is a method for class COMPARTMENT
%
%   getParamsNum(obj) return the total number of model parameters
%
%   getParamsNum(obj, mode) 
%           mode: String of type 'char': 'all', 'free', 'constant', 'dummy'
%
if nargin == 1
   mode = 'all';
end

switch mode
   case 'all'
       n = obj.nParams;

   case 'free'
       isConstantOrDummy = obj.links.getDummy | ...
                           obj.links.getConstant;
       n = sum(~isConstantOrDummy);

   case 'constant'
       isConstant = obj.links.getConstant;
       n = sum(isConstant);

   case 'dummy'
       isDummy = obj.links.getDummy();
       n = sum(isDummy);

   otherwise
       msg = strcat("Input mode '%s' is not recognised.\nAvailable ", ...
                    "options: 'all', 'free', 'constant', and 'dummy'.");
       error('MATLAB:compartment:getParamsNum', msg, mode);
end%of switch-case
end