classdef (Abstract) compartment < handle 
    % compartment 
    %
    %   compartment is an Abstract parametric model class for 
    %   data synthesis or model fitting. 
    %
    %   Note: Subclasses must implement the synthesize method. This class
    %   implements a default getJacobian method, which estimates the
    %   jacobian of the signal with respect to the model parameters
    %   numerically. For greater efficiency during model fitting,
    %   subclasses should override the numerical estimation with analytic
    %   expressions.
    %
    %   methods (public):
    %       fit             - fit model parameters to measured signals
    %       fitMultiRun     - run the optimisation using different 
    %                         initial points selected randomly.
    %       synthesize      - synthesize attenuation signals
    %       getCost         - return the Root Mean Square (RMS) error of
    %                         the fitted model
    %       getParams       - return the model parameters in a column
    %                         vector
    %       getParamsNum    - return the number of model parameters
    %       getFixedParams  - return a logical column vector stating if a
    %                         parameter is fixed (true) or not (false).
    %
    %   methods (TO BE IMPLEMENTED IN SUB-CLASSES):
    %       synthesize       - return signal from the specific model.
    %       getJacobian     - return the gradient of model wrt to its 
    %                         parameters. (analytic solution)
    %       linkFun         - define appropriate link functions for each
    %                         model parameters to enforce spicific criteria
    %                         for each parameter; ex. positivity of
    %                         diffusivity is enforced using x^2 function.
    %       invLinkFun      - define inverse link functions
    %       randomInit      - set random initial values for the parameters
    %       updateParams    - update model parameters
    %       set             - allow setting model parameters
    %
    %   properties (TO BE IMPLEMENTED IN SUB-CLASSES):
    %       fitter          - an "optimizer" object used for fitting model
    %                         parameters to data
    %       modelParams     - a column-wise numberical vector including all
    %                         model parameters
    %       fixedParams     - a column-wise boolean vector with the same
    %                         size as modelParams. if fixedParams(i) is 
    %                         false, then the i-th parmeter is included in
    %                         the optimisation.
    %       nParams         - the total number of parameters in the model
    %
    %   See also: tensor, ball, zeppelin, cylinder, ellipticalCylinder,
    %   stick, threeCompartmentModel, twoCompartmentModel
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (Constant=true, Access=protected)
        % the gyromagnetic Ratio (copied from camino source codes)
        GAMMA = 2.6751525e8; 
        
        % 60 first roots from the equation J'1(x)=0
        % J'1 is the derivative of the first order bessel function
        % (copied from camino source codes)
        Jp1ROOTS = ...
         [1.84118307861360, 5.33144196877749, 8.53631578218074, ...
          11.7060038949077, 14.8635881488839, 18.0155278304879, ...
          21.1643671187891, 24.3113254834588, 27.4570501848623, ...
          30.6019229722078, 33.7461812269726, 36.8899866873805, ...
          40.0334439409610, 43.1766274212415, 46.3195966792621, ...
          49.4623908440429, 52.6050411092602, 55.7475709551533, ...
          58.8900018651876, 62.0323477967829, 65.1746202084584, ...
          68.3168306640438, 71.4589869258787, 74.6010956133729, ...
          77.7431620631416, 80.8851921057280, 84.0271895462953, ...
          87.1691575709855, 90.3110993488875, 93.4530179063458, ...
          96.5949155953313, 99.7367932203820, 102.878653768715, ...
          106.020498619541, 109.162329055405, 112.304145672561, ...
          115.445950418834, 118.587744574512, 121.729527118091, ...
          124.871300497614, 128.013065217171, 131.154821965250, ...
          134.296570328107, 137.438311926144, 140.580047659913, ...
          143.721775748727, 146.863498476739, 150.005215971725, ...
          153.146928691331, 156.288635801966, 159.430338769213, ...
          162.572038308643, 165.713732347338, 168.855423073845, ...
          171.997111729391, 175.138794734935, 178.280475036977, ...
          181.422152668422, 184.563828222242, 187.705499575101]';
    end
    
    properties (SetAccess = 'protected')
        fitter = optimizer();           % an "optimizer" object for fitting
                                        % model parameters
    end
    
    properties (Abstract, SetAccess = 'protected')
        links;                          % vector of type LINKER that maps
                                        % constrained model parameters to
                                        % unconstrained optimisation
                                        % variables
    end
    
    properties (Abstract, Access=protected)                                
        modelParams;                    % numerical column vector including
                                        % all model parameters
        
        hyperparams;                                
        nParams;                        % total number of model parameters 
        
    end
    
    methods (Abstract, Access=protected)
        s   = synthesize(obj, scheme);
        jac = jacobian(obj, scheme);
        updateParams(obj, p);
        updateHyperparams(obj, p);
    end
    
    methods (Abstract)
        set(obj, varargin);
        rotateAxis(obj);
    end
    
    methods (Access = public)
        function s = synth(obj, scheme)
            % return synthetically generated signal
            s = obj.synthesize(scheme);
        end

        function p = fit(obj, scheme, data)
            % fti(obj, scheme, data) fit model parameters to the given
            % data.
            %

            % check consistency between data and scheme file
            if length(data)~=size(scheme, 1)
                error('MATLAB:inconsistentDimensions'            ,...
                     ['The vector d shoudl have the same number ',...
                      'of elements as the number of measurments ',...
                      'in the scheme file.\n']); 
            end
            
            % read initial values for the model parameters
            p0 = obj.modelParams;
            
            % transform model parameters to the optimisation parametrs by
            % passing through link-functions
            x0 = obj.modelToOpt(p0);
            
            % Do the optimisation using matlab built-in functions
            obj.fitter.setfeval(@(x) obj.Feval(x, scheme, data))
            obj.fitter.setfjac(@(x) obj.fjac(x, scheme));
            
            [x, cost, exitFlag] = obj.fitter.run(x0);

            p = obj.optToModel(x);
            obj.updateParams(p);
            obj.rotateAxis();
            p = [exitFlag; cost; obj.modelParams];
        end
        
        function [p, P] = fitMultiRun(obj, scheme, data, nReps)
            % fitMultiRun run the optimisation problem from different
            % starting points and return the the best-fit model parameters
            % p.
            % 
            if nargin == 3
                nReps = 50;
            end

            % check consistency between data and scheme file
            if size(data,1)~=size(scheme, 1)
                error('MATLAB:inconsistentDimensions'            ,...
                     ['The vector d shoudl have the same number ',...
                      'of elements as the number of measurments ',...
                      'in the scheme file.\n']); 
            end

            P = zeros(obj.nParams+2, nReps);
            
            obj.fitter.setfeval(@(x) obj.Feval(x, scheme, data))
            obj.fitter.setfjac(@(x) obj.fjac(x, scheme));
            parfor i = 1:nReps
                p0 = obj.randomInit();
                x0 = obj.modelToOpt(p0);
                [thisX, thisCost, exitflag] = obj.fitter.run(x0);
                thisP = obj.optToModel(thisX);
                P(:, i) = [exitflag; thisCost; thisP];
            end
            [~, iMin] = min(P(2,:));
            p = P(3:end, iMin);
            obj.updateParams(p);
            obj.rotateAxis();
            p = [P(1:2, iMin); obj.getParams()];
        end

        function cost = getCost(obj, scheme, data)
                s = obj.synthesize(scheme);
                cost = sum((data-s).^2);
        end
        
        function p = getParams(obj, option)
            if nargin == 1
                option = 'all';
            end
            
            switch option
                case 'all'
                    p = obj.modelParams;
                case 'free'
                    isDummy = obj.links.getDummy;
                    isConstant = obj.links.getConstant;
                    p = obj.modelParams(not(isDummy) & not(isConstant));
                case 'constant'
                    isConstant = obj.links.getConstant;
                    p = obj.modelParams(isConstant);
                case 'dummy'
                    isConstant = obj.links.getConstant;
                    p = obj.modelParams(isConstant);
                otherwise
                    errString = strcat("Option '%s' is not recognised.", ...
                                      "\nAvailable options are ", ...
                                      "'all', 'free', 'constant',", ...
                                      " and 'dummy'.");
                   error('MATLAB:compartment:getParams', ...
                         errString, option);
            end
        end
        
        function p = getHyperparams(obj)
            p = obj.hyperparams;
        end
        
        function links = getLinks(obj)
            links = obj.links;
        end
        
        function n = getParamsNum(obj, option)
           if nargin == 1
               option = 'all';
           end
           
           switch option
               case 'all'
                   n = length(obj.modelParams);
                   
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
                   errString = strcat("Option '%s' is not recognised.", ...
                                      "\nAvailable options are ", ...
                                      "'all', 'free', 'constant',", ...
                                      " and 'dummy'.");
                   error('MATLAB:compartment:getParamsNum', ...
                         errString, option);
           end
        end
        
        function n = getHyperparamsNum(obj)
            n = length(obj.hyperparams);
        end
        
        
        function v = getVolFraction(obj)
            v = obj.s0;
        end
        
        function addConstraint(obj, str)
            % addConstraint(obj, str) enforce various constraints on model 
            % parameters using appropriate linking functions.
            % Inputs:
            % obj:            Compartment model
            % str:            An string identifying the type of constraint 
            %                 Examples:
            %                 1) 's0 = 1' set s0 as constant with value 1
            %                 2) 'diffPar >= 1.5e-9'
            %                 3) 'diffPar <= 1.5e-9'
            %                 3) 'diffPar >= diffPerp1'
            %                 4) 'diffPerp2 = diffPerp1'
            
            % extract model parameters
            modelParamList = obj.getParamsList();
            
            % create string listing all parameters for printing if needed
            msg = [];
            for n = 1:obj.nParams-1
                msg = strcat(msg, sprintf("'%s'", modelParamList{n}), ', ');
            end
            msg =  strcat(msg, sprintf("and '%s'", modelParamList{obj.nParams}), '.');
            
            % analyse input string
            % step 1: remove all spaces
            str(isspace(str))=[];
            
            % step 2: extract parameters and replace them with pi
            [params, matches] = strsplit(str, {'+', '=', '>=', '<=', ...
                                               ')', '(', '*', '/', '>', ...
                                               '<', '^', '-', 'e-','e+'});
            for i = 1:length(params)
                thisParam = params{i};

                isNotNumeric = isnan(str2double(thisParam));
                isNotEmpty = ~isempty(thisParam);
                
                if isNotNumeric && isNotEmpty
                     % check if param1 is a valid model parameter
                    [isValidParamName, iParam] = ismember(thisParam, modelParamList);              
                    if ~isValidParamName
                        error(strcat("'%s' is not a valid parameter name.\n", ...
                               "Valid variable names for class %s are %s."),...
                               thisParam, obj.name, msg);
                    end
                    
                    params{i} = sprintf('p%d', iParam);
                end
            end
            
            str = strjoin(params, matches);
            obj.links.addConstraint(str);
                    
%             switch operator
%                 case '='
%                     if isNotNumeric
%                         thisStr = sprintf('p%d=p%d', nLink1, nLink2);
%                         obj.links.addConstraint(thisStr);
%                     else
%                         % check if numerical value is in range
%                         if value < obj.links(nLink1).lowerBound
%                             obj.links(nLink1).set('lowerBound', value);
%                             
%                         elseif value > obj.links(nLink1).upperBound
%                             obj.links(nLink1).set('upperBound', value);
%                         end
%                         thisStr = sprintf('p%d=%d', nLink1, value);
%                         obj.links.addConstraint(thisStr);
%                     end
%                     
%                 case '>='
%                     if isNotNumeric
%                         thisStr = sprintf('p%d>=p%d', nLink1, nLink2);
%                         obj.links.addConstraint(thisStr);
%                     else
%                         obj.links(nLink1).set('lowerBound', value);
%                     end
%                     
%                 case '<='
%                     if isNotNumeric
%                         thisStr = sprintf('p%d>=p%d', nLink2, nLink1);
%                         obj.links.addConstraint(thisStr);
%                     else
%                         obj.links(nLink1).set('upperBound', value);
%                     end
%                     
%                 otherwise
%                     error(strcat("'%s' is not a valid operator.\n" , ...
%                           "Valid operators are '=', '>=' or '<='."), ...
%                           operator);
%                 %\\          
%             end % of switch-case
        end %of addConstraint
        
        function constraintList = getConstraints(obj)
            constraintList = obj.links.getConstraints();
        end
        
        function p = randomInit(obj, seed)
            % randomInit initialise the model parameters randomly using 
            % plausible biophysical ranges
            if nargin == 2
                % make sure random number generation is repeatable
                rng(seed);
            end
            
            p = obj.links.init();
            obj.updateParams(p)
        end
        
        function jac = getParamsJacobian(obj, scheme)
            % return model jacobian wrt model parameters
            jac = obj.jacobian(scheme);
        end
        
        function jac = getJacobian(obj, scheme, flag)
            % getJacobian return the gradient of singal wrt to optimisation
            % parameters, i.e. model parameters passed through the linkFun.
            
            if nargin == 2
                flag = false;
            end
            
            if flag
                % compute the jacobian matrix using numerical technique.
                %\\
                EPS = 1e-7;                 % epsilon used for computation
                                            % of numerical derivatives

                p0 = obj.getParams;            % model params
                x0 = obj.links.invLink(p0);    % linking parameters
                f0 = obj.synthesize(scheme);    % signal value

                nOptParams = length(x0);
                % initialise the jac with zeros
                jac = zeros(size(scheme,1), nOptParams);
                for n = 1:nOptParams
                    x = x0; x(n) = x(n) + EPS;
                    p = obj.links.link(x);
                    obj.updateParams(p);
                    jac(:,n) = (obj.synthesize(scheme) - f0)/EPS;
                end
                obj.updateParams(p0);
                %\\
            else
                % compute jacobian matrix using analytic solution
                paramsJac = obj.jacobian(scheme);
                x = obj.links.invLink(obj.modelParams);
                linkingJac = obj.links.jacobian(x);
                jac = paramsJac * linkingJac;
            end
            %\\
        end %of getJacobian
        
        function varargout = testJacobian(obj, scheme)
            jac_num = getJacobian(obj, scheme, true);
                        
            jac_anl = getJacobian(obj, scheme, false);
            
            TOL = 1e-5;
            
            err = norm(jac_num(:)-jac_anl(:))/norm(jac_num(:));
            if err<TOL
                fprintf('Analytic solution is fine!\n');
            else
                fprintf('Analytic solution is buggy!\n');
            end
            
            if nargout == 0
                % do nothing!
            elseif nargout == 1
                varargout{1} = jac_anl;
            elseif nargout == 2
                varargout{1} = jac_anl;
                varargout{2} = jac_num;
            else
                error('Improper number of outputs.\n');
            end
        end % of testjacobian
        %\\
    end % of methods (public)
    
    methods (Access=protected)
        
        function x = modelToOpt(obj, p)
            x = obj.links.invLink(p);
        end
        
        function p = optToModel(obj, x)
            p = obj.links.link(x);
        end
        
        function f = feval(obj, x, scheme, data)
            % compute the mean squared error f=\sum (s_m - d_m)^2
            p = obj.optToModel(x);
            obj.updateParams(p);
            s = obj.synthesize(scheme);
            f = sum((data-s).^2);
        end
        
        function F = Feval(obj, x, scheme, data)
            % compute the difference between model and measurements
            % f[m] = s_m - d_m
            p = obj.optToModel(x);
            obj.updateParams(p);
            s = obj.synthesize(scheme);
            F = s - data;
        end
        
        function jac = fjac(obj, x, scheme)
            p = obj.optToModel(x);
            obj.updateParams(p);
            jac = obj.getJacobian(scheme);
        end % of fgrad
        
        function gf = fgrad(obj, x, scheme, data)
            jac = obj.fjac(x, scheme);
            F   = obj.Feval(x, scheme, data);
            gf = jac'*F*2;
        end % of fgrad
        %\\
    end % of methods (protected)
end