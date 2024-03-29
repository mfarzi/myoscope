classdef zeppelin < tensor
    % ZEPPELIN
    %
    %   A ZEPPELIN object is a prolate TENSOR where the secondary and the 
    %   tertiary diffusion eigenvalues are similar and less than the
    %   primary eigenvalue. Since the secondary (u2) and tertiary (u3)
    %   eigenvectors are free to rotate about the primary eigenvector (u1),
    %   here we assume alpha=pi/2 to address the degeneracy in the
    %   parameter space. This class only overwrites the getlinks method.
    %   This class should be used as input to MULTICOMPARTMENT class.
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
    %       Panagiotaki, E., Schneider, T., Siow, B., Hall, M.G., Lythgoe,
    %       M.F. and Alexander, D.C., "Compartment models of the diffusion
    %       MR signal in brain white matter: a taxonomy and comparison.",
    %       Neuroimage, 59(3), pp.2241-2254, 2012.
    %
    %
    %   See also: multicompartment, tensor, pancake, ball
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
        
    methods
        function obj = zeppelin()
            %ZEPPELIN Construct Function.
            %
            %   zeppelin() constructs a single-compartment model
            
            obj.name = 'zeppelin';
            
            % initialise the linking functions to map constrained model
            % parameters to unconstrained optimisation variables.
            links = obj.getLinks();
            % set diffPepr1 = diffPerp2
            links.addConstraint('diffPerp1=diffPerp2');
            % set alpha = pi/2
            links.addConstraint('alpha=pi/2');
            obj.links = links;
        end
    end % methods (public)
    %\\
end % of ZEPPELIN class
