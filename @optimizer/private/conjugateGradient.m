function [x, cost, exitFlag, output] = conjugateGradient(fval, gval, x0, options, varargin)
    % see page 278
    % Luenberger, David G., and Yinyu Ye. Linear and nonlinear programming.
    % Third Edition. Chapter 9, pp 278, Springer, 2008.
    % conjugate gradient method
    % Initialization
    x1 = x0; 
    fx1 = fval(x1, varargin{:});
    gx1 = gval(x1, varargin{:});
    Pnew = -gx1;
    
    N = length(x0); % use for restarting
    counter = 1;
    % choosing alpha1 for the line search algorithm
    deltaf = 1;
    options.alpha1 = max(2*deltaf/(gx1'*gx1),options.stopTol);
    
    %\\
    % print the result
    if options.prettyPrint
        fprintf('START of the Conjugate Gradient algorithm: \n');
        fprintf('--------------------------------------------\n');
    end

    %\\
    %--------------------------------------------------------------
    % main loop of the program starts here
    %--------------------------------------------------------------
    %\\
    % star of the loop
    for iter = 0:options.maxIter
        counter = counter + 1;
        % print the result
        if options.prettyPrint
            fprintf('iteration %2.0f --> function value = %1.8f\n',iter,fx1);
        end

        % line search algorithm   
        [alpha,flag] = wolfLineSearch(fval, gval, x1,Pnew,fx1,gx1, options, varargin{:});
        if flag == 1
            x2 = x1; break;
        end

        % make the movement
        x2 = x1 + alpha*Pnew;  
        gx2 = gval(x2, varargin{:});    
        fx2 = fval(x2, varargin{:});
        % considering stoping condition
        %if (abs(fx1-fx2) < options.stopTol)
        if norm(gx2) < options.stopTol
            break;
        else
            % find the direction
            Pold = Pnew;
            beta = (gx2'*gx2)/(gx1'*gx1); % Fletcher?Reeves method
            Pnew = -gx2 + beta*Pold; 
                       
            % resarting procedure if Pnew is not a descending direction 
            if Pnew'*gx2>=0 || counter > N
                Pnew = -gx2;
                counter = 1;
            end
            
            deltaf = fx1 - fx2;
            options.alpha1 = max(2*deltaf/(gx2'*gx2),options.stopTol);
            
            % updating the old x1 with new x2
            fx1 = fx2; gx1 = gx2;
            x1 = x2;
        end
    end % of for loop
    %
    % print the result
    if options.prettyPrint
        fprintf('END of the algorithm: \n');
        fprintf('--------------------------------------------\n');
    end
    x = x2;
    cost = fval(x, varargin{:});
    output.nIter = iter;
    exitFlag = 1;
end % of conjugate gradient method