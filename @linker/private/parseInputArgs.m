function args = parseInputArgs(varargin)
    % parseInputArgs is a private method for class LINKER
     parser = inputParser;
     parser.CaseSensitive = false;
     parser.addRequired('name', @(v) ischar(v));
     parser.addOptional('type', 'unbounded', @(v) ischar(v) && ...
         (strcmp(v, 'unbounded')||strcmp(v, 'bounded')));
     parser.addOptional('minimum', [], @(v) isnumeric(v) && ...
         isscalar(v));
     parser.addOptional('maximum', [], @(v) isnumeric(v) && ...
         isscalar(v));
     parser.addParameter('initRange', [-1 1], ...
         @(v) validateattributes(v, {'numeric'}, {'size', [1,2]})); 
     
     parser.parse(varargin{:});
     args = parser.Results;
     
      % check maximum and minimums are set properly
     if strcmp(args.type, 'bounded')
         assert(~isempty(args.minimum), 'MATLAB:linker',...
             'For a bounded link type, minimum (3rd argument) cannot be empty.');
         assert(~isempty(args.maximum), 'MATLAB:linker',...
             'For a bounded link type, maximum (4th argument) cannot be empty.');
         assert(args.maximum>args.minimum, 'MATLAB:linkder',...
             ['For a bounded link type, maximum (4th argument)',...
             ' must be greater than the minimum (3rd argument).']);         
     end
     
     if strcmp(args.type, 'unbounded')
         assert(isempty(args.minimum), 'MATLAB:linker',...
             'For an unbounded link type, minimum (3rd argument) must be empty.');
         assert(isempty(args.maximum), 'MATLAB:linker',...
             'For an unbounded link type, maximum (4th argument) must be empty.');
     end
end