classdef ellipticalCylinder < compartment
    % ELLIPTICALCYLINDER 
    % (Elliptical Cylinder Model with Gaussian Phase Distribution)
    %
    %   An ELLIPTICALCYLINDER object is a basic COMPARTMENT object 
    %   represenitng particles diffusion in an elliptical cylinder. This
    %   class encapsulates the model parameters and basic methods for 
    %   fitting the parameters to diffusion weigthed MR signals or 
    %   synthesize signals for a given diffusion scheme.
    %
    %   For mathematical background see
    %       Vangelderen, P., DesPres, D., Vanzijl, P.C.M. and Moonen, C.,
    %       "Evaluation of restricted diffusion in cylinders.
    %       Phosphocreatine in rabbit leg muscle.", Journal of Magnetic
    %       Resonance, Series B, 103(3), pp.255-260, 1994.
    %
    %   properties:
    %       name              - class name
    %       s0                - b0-signal with no diffusion weight or
    %                           signal attenuation at b-val=0 
    %       diff              - diffusivity along the cylinder axis [m^2/s]
    %       r1                - radius of ellipsoid along the major axis
    %                           [m]
    %       r2                - radius of ellipsoid along the minor axis
    %                           [m]
    %       theta             - elevation angle; the angle between the
    %                           cylinder axis (n1) and the z-axis. [radian]
    %       phi               - azimuth angle; the angle between the x-axis
    %                           and the n1 projection onto the xy-plane. 
    %                           [radian]
    %       alpha             - the angle between the major ellipsoid axis
    %                           and the n1 projection onto the cross 
    %                           sectional elliposid plain [radian]
    %       fitter            - an "optimizer" object to fit model
    %                           parameters
    %
    %    methods:
    %       fit               - fit model parameters to diffusion signals.
    %                           see the class COMPARTMENT.
    %       fitMultiRun       - fit model parameters using different set of
    %                           start points. see the class COMPARTMENT.
    %       synthesize        - synthesize signals for a diffusion scheme
    %                           see the class COMPARTMENT for definition.
    %       randomInit        - randomly initialise model parameters.
    %       set               - allow setting model parameters.
    %       fixParams         - Set specific parameters as constant during
    %                           model fitting. By default, all parameters
    %                           are set to be variables.
    %       getCost           - return the Root Mean Square (RMS) error of
    %                           the fitted model. see class COMPARTMENT.
    %       getParams         - return the model parameters in a column
    %                           vector. see class COMPARTMENT.
    %       getParamsNum      - return the number of model parameters.
    %                           see class COMPARTMENT.
    %       getFixedParams    - return a logical column vector stating if a
    %                           parameter is fixed (true) or not (false).
    %       getCylinderAxis   - return the normal vector (n1) giving the
    %                           cylinder axis orientation.
    %
    %   See also: compartment, cylinder, stick
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess = 'protected')
        name  = 'ellipticalCylinder'       % class name
        s0    = 1;                  % b0-signal with no diffusion weight
        
        diffPar  = 1e-9;            % diffusivity along the cylinder axis
                                    % [m^2/s]
                                       
        r1    = 10e-6;              % Radius of ellipsoid along the major
                                    % axis [m]
                                       
        r2    = 5e-6;               % Radius of ellipsoid along the minor
                                    % axis [m]
                                       
        theta = 0;                  % elevation angle; the angle between 
                                    % the cylinder axis (n1) and the z-axis 
                                    % [radian]
                                       
        phi   = 0;                  % azimuth angle; the angle between the
                                    % x-axis and the n1 projection  onto
                                    % the xy-plane. [radian]

        alpha = pi/2;               % the angle between the major ellipsoid
                                    % axis and the n1 projection onto the
                                    % cross sectional elliposid plain 
                                    % [radian]     
                                    
        links     = [];             % vector of type LINKER that maps 
                                    % constrained model parameters to
                                    % unconstrained optimisation
                                    % variables                            
    end
    
    
    properties (Access=protected)
        nParams = 7;                % number of model parameters
        modelParams  = [];          % model parameters 
                                    % [s0; diffPar; r1; r2; theta; phi; alpha]   
        hyperparams = [];                            
    end
    
    methods 
        function obj = ellipticalCylinder(varargin) 
            %ELLIPTICALCYLINDERGPD Construct Function.
            %   ellipticalCylinder() construct an object with default 
            %   sets of parameters.
            %
            %   ellipticalCylinder(params) construct an object with 
            %   given initial params; [diffPar; r1; r2; theta; phi; alpha]
            %
            %   ellipticalCylinder(..., <Parameter>, <Value>, ...) allow
            %   passing model parameters as parameter-value pairs. If
            %   params is provided as first argument, the param-value pairs
            %   will be ignored.
            % 
            obj.modelParams = [obj.s0; obj.diffPar; obj.r1; obj.r2; ...
                               obj.theta; obj.phi; obj.alpha];               
            
            % define linkers for the zeppelin class
            obj.links = [linker('type', 'cos'   , 'lowerBound', 0, 'upperBound', 1, 'initRange', [0, 1]);  ...
                         linker('type', 'cos'   , 'lowerBound', 0.1e-9, 'upperBound', 3e-9, 'initRange', [1e-10, 3e-9]);  ...
                         linker('type', 'cos'   , 'lowerBound', 0.1e-6, 'upperBound', 20e-6, 'initRange', [0.1e-6, 20e-6]);  ...
                         linker('type', 'cos'   , 'lowerBound', 0.1e-6, 'upperBound', 20e-6, 'initRange', [0.1e-6, 20e-6]);  ...
                         linker('type', 'linear', 'initRange', [0, pi/2]);  ...
                         linker('type', 'linear', 'initRange', [0, 2*pi]);
                         linker('type', 'linear', 'initRange', [0, pi])];               
            obj.links.setName(obj.getParamsList);
            
            if nargin >0
                if ~isa(varargin{1}, 'char')
                    obj.set('params', varargin{:});
                    if ismember('params', varargin{:})
                        warning(['The first optional input is ignored', ...
                                 ' since model parameters are set', ...
                                 ' using value pair for "params".\n']);
                    end
                else
                    obj.set(varargin{:});
                end
            end
        end % of constructor    
        
        function set(obj, varargin)
            % set(obj, varargin) allow changing model parameters for an
            % existing object.
            %
            % obj.set() do nothing!
            %
            % obj.set(..., <Parameter>, <Value>, ...) 
            % allow passing model parameters as parameter-value pairs. 
            %
            % NOTE: If 'params' are provided, the param-value pairs will be
            % ignored.
            %
            if nargin==1
                % >> obj = obj.set()
                % Do Nothing!
            else
                % Parse the parameter-value pairs
                p = inputParser;
                p.CaseSensitive = false;
                p.addParameter('params'        , obj.modelParams   , @obj.validationFCN_params);
                p.addParameter('s0'            , obj.s0            , @obj.validationFCN_s0);                        
                p.addParameter('diffPar'       , obj.diffPar       , @obj.validationFCN_diffPar);
                p.addParameter('r1'            , obj.r1            , @obj.validationFCN_r1);
                p.addParameter('r2'            , obj.r2            , @obj.validationFCN_r2);
                p.addParameter('theta'         , obj.theta         , @obj.validationFCN_theta);
                p.addParameter('phi'           , obj.phi           , @obj.validationFCN_phi);
                p.addParameter('alpha'         , obj.alpha         , @obj.validationFCN_alpha);   
                p.addParameter('fitter'        , obj.fitter        , @obj.validationFCN_fitter);
                
                p.parse(varargin{:});              
                
                % update the fitter
                obj.fitter = p.Results.fitter;

                paramsIsProvided = not(ismember('params', p.UsingDefaults));
                
                if paramsIsProvided
                    obj.s0        = p.Results.params(1);
                    obj.diffPar   = p.Results.params(2);
                    obj.r1        = p.Results.params(3);
                    obj.r2        = p.Results.params(4);
                    obj.theta     = p.Results.params(5);
                    obj.phi       = p.Results.params(6);
                    obj.alpha     = p.Results.params(7);
                    
                    % warning if parameters are overwritten
                    paramList = {'s0', 'diffPar', 'r1', 'r2', ...
                                 'theta', 'phi', 'alpha'};
                    id = find(~(ismember(paramList, p.UsingDefaults)));
                    if length(id) == 1
                        warning(['Parameters are already provided using ', ...
                                 'vector format. The input for "%s" ', ...
                                 'is ignored.'], paramList{id});
                    elseif length(id) > 1
                        str = sprintf('"%s"', paramList{id(1)});
                        for i = 2:length(id)
                            if i == length(id)
                                str = strcat(str, sprintf(' and "%s"', paramList{id(i)}));
                            else
                                str = strcat(str, sprintf(', "%s"', paramList{id(i)}));
                            end
                        end
                        warning(['Parameters are already provided using ', ...
                                 'vector format. The inputs for %s ', ...
                                 'are ignored.'], str);
                    end
                %\\    
                else % use param-value pairs to set the parameters
                    obj.s0        = p.Results.s0;
                    obj.diffPar   = p.Results.diffPar;
                    obj.r1        = p.Results.r1;
                    obj.r2        = p.Results.r2;
                    obj.theta     = p.Results.theta;
                    obj.phi       = p.Results.phi;
                    obj.alpha     = p.Results.alpha;
                end
        
                obj.modelParams = [obj.s0; obj.diffPar; obj.r1; obj.r2; ...
                                   obj.theta; obj.phi; obj.alpha];             
            end % of if nargin==0
        end % of set
        
        function rotateAxis(obj)
           % rotateAxis remove the ambiguity in estimation of parameters
           % theta, phi, and alpha by rotating the cooriante system of the
           % elliptical cylinder appropriately.
           % 
           % The orthonormal coordiante system V = [v1, v2, v3] should 
           % rotate such that v1 is parallel to the cylinder axis, v2 is 
           % parallel to the major diameter of the ellipsoid, and v3 is 
           % parallel to the minor diameter of the elliposid.
           %
           % More over, the direction of v1 or -v1 is selected such that
           % the angle theta is always in range [0, pi/2].
           % 
           % Finally, the direction v2 or -v2 is selected such that the
           % angle alpha is always between [0, phi].
           %
           % see also: getUnitFrame, getEulerAngles
            
           V = getUnitFrame(obj.theta, obj.phi, obj.alpha);
            
           if obj.r1<obj.r2
               V = V(:, [1,3,2]);
               V(:,3) = -V(:,3);
               r = obj.r1;
               obj.r1 = obj.r2;
               obj.r2 = r;
               warning('MATLAB:ellipticalCylinder:rotateAxis', ...
                       ' Constraint on r1 and r2 is not set properly.');
           end

           [obj.theta, obj.phi, obj.alpha] = getEulerAngles(V);
            
           % make sure theta < pi/2
           if obj.theta>pi/2
               obj.theta = pi - obj.theta;
               obj.phi = pi + obj.phi;
               obj.alpha = -obj.alpha;
           end

           obj.phi = mod(obj.phi, 2*pi);
            
           % make sure alpha is in rage [0, pi]
           obj.alpha = mod(obj.alpha, pi);
           
           obj.modelParams(2:end) = [obj.diffPar; obj.r1;...
                                     obj.r2; obj.theta; ...
                                     obj.phi; obj.alpha];
        end % of rotateAxis
        
        function n = getCylinderAxis(obj)
            % compute the normal vector parallel to the cylinder axis
            n = [cos(obj.phi)*sin(obj.theta); sin(obj.phi)*sin(obj.theta); cos(obj.theta)];
        end
        
        function N = getAxes(obj)
            n1 = [cos(obj.phi)*sin(obj.theta); sin(obj.phi)*sin(obj.theta); cos(obj.theta)];
            
            rst = sin(obj.theta+pi/2);
            n2_zrot = [  rst * cos(obj.phi); ...
                         rst * sin(obj.phi); ...
                        cos(obj.theta+pi/2)];
        
            % rotate around n1 by alpha
            n2 = rodrigues_rot(n2_zrot, n1, obj.alpha);
            n3 = cross(n1, n2);
            N = [n1, n2, n3];
        end
    end % of methods (public)
        
    methods (Access = protected)
        %\\
        function s = synthesize(obj, scheme)
            % synthesize(obj, scheme) synthesize signals for diffusion
            % particles in an elliptical cylinder. The signal is modeled as
            % the product of the signals parallel and perpendicular to the 
            % cylinder axis
            % 
            
            G = [scheme.x, scheme.y, scheme.z].*scheme.G_mag;
            N = obj.getAxes(); n1 = N(:,1); n2 = N(:,2); n3 = N(:,3);
            
            % comput the gradient along and perpendicular
            % to the cylinder axis 
            Gpar_mag2 = (G*n1).^2; 
            Gperp1_mag2 = (G*n2).^2;
            Gperp2_mag2 = (G*n3).^2;
            
            % compute bPar
            bPar = (scheme.DELTA-scheme.delta/3).*(obj.GAMMA*scheme.delta).^2*obj.diffPar;
            
            % compute bPerp1
            beta = obj.Jp1ROOTS/obj.r1;          
            DB2 = obj.diffPar*beta.^2;
            nom = 2*scheme.delta*DB2' - 2  ...
                + 2*exp(-scheme.delta*DB2') ...
                + 2*exp(-scheme.DELTA*DB2') ...
                - exp(-(scheme.DELTA-scheme.delta)*DB2') ...
                - exp(-(scheme.DELTA+scheme.delta)*DB2');

            denom = obj.diffPar^2*beta.^6.*(obj.Jp1ROOTS.^2-1);
            
            bPerp1 = 2*obj.GAMMA^2*nom*(1./denom);
            
            % compute bPerp2
            beta = obj.Jp1ROOTS/obj.r2;          
            DB2 = obj.diffPar*beta.^2;
            nom = 2*scheme.delta*DB2' - 2  ...
                + 2*exp(-scheme.delta*DB2') ...
                + 2*exp(-scheme.DELTA*DB2') ...
                - exp(-(scheme.DELTA-scheme.delta)*DB2') ...
                - exp(-(scheme.DELTA+scheme.delta)*DB2');

            denom = obj.diffPar^2*beta.^6.*(obj.Jp1ROOTS.^2-1);
            
            bPerp2 = 2*obj.GAMMA^2*nom*(1./denom);
            
            
            % compute signal as the product of the value along and
            % perpendicular to the cylinder axis
            s = obj.s0 * exp(-bPar  .*Gpar_mag2)   ...
                      .* exp(-bPerp1.*Gperp1_mag2) ...
                      .* exp(-bPerp2.*Gperp2_mag2);
        end % of synthesize
        
        function jac = jacobian(obj, scheme)
            % getJacobian return the gradient of singal wrt to optimisation
            % parameters, i.e. model parameters passed through the linkFun.
            %
            % NOTE: IF YOU CHANGE linkFun or invLinkFun, MAKE SURE TO EDIT
            % THIS CODE APPROPRIATELY AS WELL.
            %
            
            
            nScheme = size(scheme, 1);
            
            % initialise the jac with zeros
            jac = zeros(nScheme, obj.nParams);
            
            N = obj.getAxes();                  
            n1 = N(:,1);                        % cylinder axis
            n2 = N(:,2);                        % perpendicular axis 1
            n3 = N(:,3);                        % perpendicular axis 2
           
            G = [scheme.x, scheme.y, scheme.z].*scheme.G_mag; % gradient vector
            Gpar_mag2 = (G*n1).^2; 
            Gperp1_mag2 = (G*n2).^2;
            Gperp2_mag2 = (G*n3).^2;
            
            % compute bPar
            bPar = (scheme.DELTA-scheme.delta/3).*(obj.GAMMA*scheme.delta).^2*obj.diffPar;
            
            % compute bPerp1
            beta_r1 = obj.Jp1ROOTS/obj.r1;          
            B2_r1 = beta_r1.^2;
            
            nom1 = 2*obj.diffPar*scheme.delta*B2_r1' - 2   ...
                 + 2*exp(-obj.diffPar*scheme.delta*B2_r1') ...
                 + 2*exp(-obj.diffPar*scheme.DELTA*B2_r1') ...
                 - exp(-obj.diffPar*(scheme.DELTA-scheme.delta)*B2_r1') ...
                 - exp(-obj.diffPar*(scheme.DELTA+scheme.delta)*B2_r1');

            denom1 = obj.diffPar^2*beta_r1.^6.*(obj.Jp1ROOTS.^2-1);
            
            bPerp1 = 2*obj.GAMMA^2*nom1*(1./denom1);

            % compute bPerp2
            beta_r2 = obj.Jp1ROOTS/obj.r2;          
            B2_r2 = beta_r2.^2;

            nom2 = 2*obj.diffPar*scheme.delta*B2_r2' - 2  ...
                 + 2*exp(-obj.diffPar*scheme.delta*B2_r2') ...
                 + 2*exp(-obj.diffPar*scheme.DELTA*B2_r2') ...
                 -   exp(-obj.diffPar*(scheme.DELTA-scheme.delta)*B2_r2') ...
                 -   exp(-obj.diffPar*(scheme.DELTA+scheme.delta)*B2_r2');

            denom2 = obj.diffPar^2*beta_r2.^6.*(obj.Jp1ROOTS.^2-1);
            
            bPerp2 = 2*obj.GAMMA^2*nom2*(1./denom2);
            
            % compute feval with the current parameters
            f0 = obj.s0*exp(-bPar.*Gpar_mag2).*exp(-bPerp1.*Gperp1_mag2)...
                      .*exp(-bPerp2.*Gperp2_mag2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%       Gradient wrt to s0-jac(:,1)        %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % gradient wrt to s0
            jac(:,1) = f0/obj.s0;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%       Gradient wrt to diffPar-jac(:,2)      %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % gbPar_gDiffPar
            gbPar_gDiffPar = bPar/obj.diffPar;
            
            % gbPerp1_gDiffPar
            gNom1_gDiffPar = 2*scheme.delta*B2_r1' ...
                        - 2*exp(-obj.diffPar*scheme.delta*B2_r1').*(scheme.delta*B2_r1') ...
                        - 2*exp(-obj.diffPar*scheme.DELTA*B2_r1').*(scheme.DELTA*B2_r1') ...
                        +   exp(-obj.diffPar*(scheme.DELTA-scheme.delta)*B2_r1').*((scheme.DELTA-scheme.delta)*B2_r1')...
                        +   exp(-obj.diffPar*(scheme.DELTA+scheme.delta)*B2_r1').*((scheme.DELTA+scheme.delta)*B2_r1');
            
            gDenom1_gDiffPar = 2*obj.diffPar*beta_r1.^6.*(obj.Jp1ROOTS.^2-1);
            tmp = bsxfun(@times, gNom1_gDiffPar, denom1') - bsxfun(@times, nom1, gDenom1_gDiffPar');
            tmp = bsxfun(@rdivide, tmp, denom1'.^2);
            gbPerp1_gDiffPar = sum(tmp, 2)*2*obj.GAMMA^2;
     
            % gbPerp2_gDiffPar
            gNom2_gDiffPar = 2*scheme.delta*B2_r2' ...
                        - 2*exp(-obj.diffPar*scheme.delta*B2_r2').*(scheme.delta*B2_r2') ...
                        - 2*exp(-obj.diffPar*scheme.DELTA*B2_r2').*(scheme.DELTA*B2_r2') ...
                        +   exp(-obj.diffPar*(scheme.DELTA-scheme.delta)*B2_r2').*((scheme.DELTA-scheme.delta)*B2_r2')...
                        +   exp(-obj.diffPar*(scheme.DELTA+scheme.delta)*B2_r2').*((scheme.DELTA+scheme.delta)*B2_r2');
            
            gDenom2_gDiffPar = 2*obj.diffPar*beta_r2.^6.*(obj.Jp1ROOTS.^2-1);
            tmp = bsxfun(@times, gNom2_gDiffPar, denom2') - bsxfun(@times, nom2, gDenom2_gDiffPar');
            tmp = bsxfun(@rdivide, tmp, denom2'.^2);
            gbPerp2_gDiffPar = sum(tmp, 2)*2*obj.GAMMA^2;
            
            % gradient wrt diffPar       
            gf_gDiffPar = -(gbPar_gDiffPar.*Gpar_mag2     + ...
                         gbPerp1_gDiffPar.*Gperp1_mag2 + ...
                         gbPerp2_gDiffPar.*Gperp2_mag2).*f0;
            jac(:,2) = gf_gDiffPar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%       Gradient wrt to r1-jac(:,3)        %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % gbPar_gR1 = 0;
            
            % gbPerp1_gR1
            gNom1_gBm = 4*obj.diffPar*scheme.delta*beta_r1' ...
                      - 2*exp(-obj.diffPar*scheme.delta*B2_r1').*(2*obj.diffPar*scheme.delta*beta_r1') ...
                      - 2*exp(-obj.diffPar*scheme.DELTA*B2_r1').*(2*obj.diffPar*scheme.DELTA*beta_r1') ...
                      +   exp(-obj.diffPar*(scheme.DELTA-scheme.delta)*B2_r1').*(2*obj.diffPar*(scheme.DELTA-scheme.delta)*beta_r1') ...
                      +   exp(-obj.diffPar*(scheme.DELTA+scheme.delta)*B2_r1').*(2*obj.diffPar*(scheme.DELTA+scheme.delta)*beta_r1');
            
            gDenom1_gBm = 6*obj.diffPar^2*beta_r1.^5.*(obj.Jp1ROOTS.^2-1);
            
            tmp = bsxfun(@times, gNom1_gBm, denom1') - bsxfun(@times, nom1, gDenom1_gBm');
            gbPerp1_gBm = 2*obj.GAMMA^2*bsxfun(@rdivide, tmp, denom1'.^2);
            
            % gbPerp1_gR1
            gBm_gR1 = -obj.Jp1ROOTS/obj.r1^2;
            gbPerp1_gR1 = gbPerp1_gBm * gBm_gR1;
            
            % gradient wrt R1
            gf_gR1 = -(gbPerp1_gR1.*Gperp1_mag2).*f0;
            jac(:,3) = gf_gR1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%       Gradient wrt to r2-jac(:,4)        %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % gbPerp2_gR
            gNom2_gBm = 4*obj.diffPar*scheme.delta*beta_r2' ...
                      - 2*exp(-obj.diffPar*scheme.delta*B2_r2').*(2*obj.diffPar*scheme.delta*beta_r2') ...
                      - 2*exp(-obj.diffPar*scheme.DELTA*B2_r2').*(2*obj.diffPar*scheme.DELTA*beta_r2') ...
                      +   exp(-obj.diffPar*(scheme.DELTA-scheme.delta)*B2_r2').*(2*obj.diffPar*(scheme.DELTA-scheme.delta)*beta_r2') ...
                      +   exp(-obj.diffPar*(scheme.DELTA+scheme.delta)*B2_r2').*(2*obj.diffPar*(scheme.DELTA+scheme.delta)*beta_r2');
            
            gDenom2_gBm = 6*obj.diffPar^2*beta_r2.^5.*(obj.Jp1ROOTS.^2-1);
            
            tmp = bsxfun(@times, gNom2_gBm, denom2') - bsxfun(@times, nom2, gDenom2_gBm');
            gbPerp2_gBm = 2*obj.GAMMA^2*bsxfun(@rdivide, tmp, denom2'.^2);
            
            % gbPerp2_gR2
            gBm_gR2 = -obj.Jp1ROOTS/obj.r2^2;
            gbPerp2_gR2 = gbPerp2_gBm * gBm_gR2;
            
            % gradient wrt R2
            gf_gR2 = -(gbPerp2_gR2.*Gperp2_mag2).*f0;
            jac(:,4) = gf_gR2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%      Gradient wrt to theta-jac(:,5)      %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gf_gN1 = repmat(-2*bPar.*(G*n1).*f0,1,3).*G;
            gf_gN2 = repmat(-2*bPerp1.*(G*n2).*f0,1,3).*G;
            gf_gN3 = repmat(-2*bPerp2.*(G*n3).*f0,1,3).*G;
            
            gN1_gTheta = [ cos(obj.phi)*cos(obj.theta); ...
                           sin(obj.phi)*cos(obj.theta); ...
                          -sin(obj.theta)];
                      
            gN2_gTheta = [-cos(obj.alpha)*cos(obj.phi)*sin(obj.theta); ...
                          -cos(obj.alpha)*sin(obj.phi)*sin(obj.theta); ...
                          -cos(obj.alpha)*cos(obj.theta)];
                     
            gN3_gTheta = [ sin(obj.theta)*cos(obj.phi)*sin(obj.alpha);...
                           sin(obj.theta)*sin(obj.phi)*sin(obj.alpha);...
                           cos(obj.theta)*sin(obj.alpha)];
                       
            jac(:,5) = gf_gN1*gN1_gTheta + gf_gN2*gN2_gTheta + gf_gN3*gN3_gTheta;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%      Gradient wrt to phi-jac(:,6)         %%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gN1_gPhi   = [-sin(obj.phi)*sin(obj.theta); ...
                           cos(obj.phi)*sin(obj.theta); ...
                           0];
                       
            gN2_gPhi   = [-sin(obj.phi)*cos(obj.alpha)*cos(obj.theta) - cos(obj.phi)*sin(obj.alpha); ...          
                           cos(obj.phi)*cos(obj.alpha)*cos(obj.theta) - sin(obj.phi)*sin(obj.alpha); ...
                           0];
                       
            gN3_gPhi   = [ sin(obj.phi)*cos(obj.theta)*sin(obj.alpha) - cos(obj.phi)*cos(obj.alpha);...
                          -cos(obj.phi)*cos(obj.theta)*sin(obj.alpha) - sin(obj.phi)*cos(obj.alpha);...
                         0];    
                     
            jac(:,6) = gf_gN1*gN1_gPhi   + gf_gN2*gN2_gPhi   + gf_gN3*gN3_gPhi;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%      Gradient wrt to alpha-jac(:,7)       %%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gN2_gAlpha = [-sin(obj.alpha)*cos(obj.theta)*cos(obj.phi) - cos(obj.alpha)*sin(obj.phi);...
                          -sin(obj.alpha)*cos(obj.theta)*sin(obj.phi) + cos(obj.alpha)*cos(obj.phi);
                           sin(obj.alpha)*sin(obj.theta)];                     
           
            gN3_gAlpha = [-cos(obj.alpha)*cos(obj.theta)*cos(obj.phi)+sin(obj.alpha)*sin(obj.phi);...
                          -cos(obj.alpha)*cos(obj.theta)*sin(obj.phi)-sin(obj.alpha)*cos(obj.phi);...
                           cos(obj.alpha)*sin(obj.theta)];
           
            jac(:,7) = gf_gN2*gN2_gAlpha + gf_gN3*gN3_gAlpha;                 
        end % of jacobian
        
        function updateParams(obj, p)
            obj.s0 = p(1);
            obj.diffPar = p(2);  % diffusivity along the cylinder axis [s/m^2]
            obj.r1 = p(3);    % Radius of the cylinder [m]
            obj.r2 = p(4);    % Radius of the cyliner [m]
            obj.theta = p(5); % the angle between cylinder axis n and the z-axis
            obj.phi = p(6);   % the angle between projection of n onto the xy-plan
                              % and the x-axis 
            obj.alpha = p(7);
            obj.modelParams = p;
        end    
        
        function updateHyperparams(obj, p)
            % do nothing
        end
    end % of methods (protected)
    
    methods (Static)
        function paramList = getParamsList()
            paramList = {'s0', 'diffPar', 'r1', 'r2', 'theta', 'phi', 'alpha'};
        end
        
        function hyperparamsList = getHyperparamsList()
            hyperparamsList = {};
        end
    end
    
    methods (Access = 'private')
        function validationFCN_params(obj, v)
            assert(numel(v) == obj.nParams && isnumeric(v),...
                   sprintf(['Value must be a numeric column or row',...
                            ' vector of size %d.\n'], obj.nParams));
                        
            assert(v(1) > 0, 'Value must be positive for "s0".');
            assert(v(2) > 0, 'Value must be positive for "diffPar".');
            assert(v(3) > 0, 'Value must be positive for "r1".');
            assert(v(4) > 0, 'Value must be positive for "r2".');
        end
        
        function validationFCN_s0(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_diffPar(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_r1(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_r2(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_theta(obj, v)
            assert(isnumeric(v) && isscalar(v), ...
                   'Value must be scalar and numeric.');
        end
        
        function validationFCN_phi(obj, v)
            assert(isnumeric(v) && isscalar(v), ...
                   'Value must be scalar and numeric.');
        end
        
        function validationFCN_alpha(obj, v)
            assert(isnumeric(v) && isscalar(v), ...
                   'Value must be scalar and numeric.');
        end
        
        function validationFCN_fitter(obj, v)
            assert(isa(v, 'optimizer'), ...
                   'Value must be of type "optimizer".');
        end
        
        function validationFCN_fixParamsVec(obj, v)
            assert(numel(v) == obj.nParams && islogical(v),...
                   sprintf(['Value must be a logical column or row',...
                            ' vector of size %d.\n'], obj.nParams));
        end
        
        function validationFCN_fixParams(obj, v)
            assert(isscalar(v) && islogical(v),...
                   'Value must be scalar and logical.');
        end
    end % of methods (private)
end % of ellipticalCylinder class
