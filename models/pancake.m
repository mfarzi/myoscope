classdef pancake < tensor
    % PANCAKE
    %
    %   A PANCAKE object is an oblate TENSOR where the primary and the 
    %   secondary diffusion eigenvalues are similar and greater than the
    %   tertiary eigenvalue. Since the primary (u1) and secondary (u2)
    %   eigenvectors are free to rotate about the third eigenvector (u3),
    %   here we assume alpha=pi/2 to address the degeneracy in the
    %   parameter space. This class only overwrites the getlinks method.
    %
    %   Model Parameters:
    %       s0                - Normalised b0-signal 
    %       diffPar           - parallel diffusivity along the primary axis
    %                           [m^2/s]
    %       diffPerp1         - perpendicular diffusivity along the
    %                           sceondary axis [m^2/s]
    %       diffPerp2         - perpendicualr diffusivity along the
    %                           tertiary axis [m^2/s]
    %       theta             - elevation angle; the angle between the
    %                           first eigen vector (v1) and the z-axis.
    %                           [radian]
    %       phi               - azimuth angle; the angle between the x-axis
    %                           and the v1 projection onto the xy-plane. 
    %                           [radian]
    %       alpha             - the angle between the second eigen vector
    %                           and the v1 rotated by pi/2 around the 
    %                           z-axis. [radian]
    %
    %   For mathematical background see
    %       Basser, P.J., Mattiello, J. and LeBihan, D., "MR diffusion 
    %       tensor spectroscopy and imaging.", Biophysical journal, 66(1),
    %       pp.259-267, 1994.
    %
    %
    %   See also: tensor, ball, zeppelin, and multicompartment
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk   
    
    methods
        function obj = pancake()
            %PANCAKE Construct Function.
            %
            %   pancake() constructs a single-compartment model
            
            obj.name = 'pancake';
            
            % initialise the linking functions to map constrained model
            % parameters to unconstrained optimisation variables.
            links = obj.getLinks();
            % set diffPar = diffPerp1
            links.addConstraint('diffPar=diffPerp1');
            % set alpha = pi/2
            links.addConstraint('alpha=pi/2');
            obj.links = links;
        end
    end % methods (public)
end% of PANCAKE class
