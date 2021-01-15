classdef scheme < handle
    % SCHEME
    %
    %   A SCHEME object facilitate the management of different diffusion
    %   scheme, i.e. stejskal-tanner or b-vecotr. 
    %
    %   properties:
    %       directions        - Diffusion directionsmethod for linking function
    %       upperBound        - upper bound for parameter p
    %       lowerBound        - lower bound for parameter p
    %       x                 - memory to store current x
    %       p                 - memory to store current p
    %
    %   methods (public):
    %       link              - forward transform p = F(x)  
    %       invLink           - backward transform x = F^(-1)(p)
    %       grad              - gradient of forward transform wrt to x
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess = 'private')
        type = 'STEJSKAL-TANNER';      
    end
    
    methods 
        function obj = scheme(varargin)
            %SCHEME Construct Function.
            %
            %   scheme(bw) construct an empty object
        end
        write(obj, filename);
    end
    
    methods (Static)
        [obj, gVec] = read(filename);
    end
end%of class scheme
    