 function [alpha_k,flag] = wolfLineSearch(fval, gval, x_k, P_k, f0, gf0, options, varargin)
    % line search algorithms based on strong Wolfe conditions
    % Performs search for point satisfying the two strong Wolfe conditions:
    %
    %	f(x_k + alpha*P_k) <= f(x_k) + c1*alpha* P_k'*gf(x_k) (WC1)
    %	|P_k'.gf(x_k + alpha*P_k)| <= c2|P_k'*gf(x_k)|		  (WC2)
    %
    % Input:
    %       f         : handle function
    %       gf        : handle function for the gradient of f
    %       x_k       : current solution at iteration k for the main loop
    %       P_k       : direction for minimizing function f
    %       f0        :      f(x_k)
    %       gf0       :     gf(x_k)
    %       OptPrm    :    [struct]
    %              .alphaMax --> wolf line search: parameter used in bracketing
    %                            phase to estimate a_j and b_j at [0,alpha_max]
    %              .c1 --> parameter for WC1 with typical value of 0.01 or 1e-4
    %              .c2 --> parameter for WC2 with typical value of 0.8 or 0.9
    %              .alpha1 --> first guess for start of the algorithm
    %
    % varargin:  other parameters that may be used by function handle f and gf
    %
    % Output:
    %       alpha_k  : parameter alpha for iteration k of the main optimization
    %                  loop
    %
    %       flag     : flag indicates that whether P_k is not a good direction
    %                  or decreasing function f is not further possible in
    %                  this direction
    % 
    % Mohsen Farzi
    % Email: mhnfarzi@gamil.com
    %==========================================================================
    %==========================================================================
    alpha_max = options.alphaMax;
    alpha1    = options.alpha1;
    c1        = options.c1;
    c2        = options.c2;

    % initializing
    flag = 0;
    phi_prime_0 = gf0'*P_k;
    phi_0 = f0;
    tau1 = 10;

    % check whether P_k is a descent direction
    if phi_prime_0 >= 0
        fprintf('Error: Need a descent direction\n');
        flag = 1;
        alpha_k = 0;
        return;
    end

    alpha_new = alpha1;
    alpha_old = 0;
    phi_alpha_old = phi_0;
    phi_prime_alpha_old = phi_prime_0;
    % bracketing phase
    while(1)
        % Evaluate phi(alpha_i)
        %[phi_alpha_new, flag] = fval(x_k + alpha_new*P_k, varargin{:});
        phi_alpha_new = fval(x_k + alpha_new*P_k, varargin{:});

        if isnan(phi_alpha_new)
            % function vaue is nan
            % stop the optimisation here!
            alpha_k = 0;
            return;
        end

        if (phi_alpha_new > phi_0 + c1*alpha_new*phi_prime_0) || (phi_alpha_new >= phi_alpha_old)
            a = alpha_old;
            phi_a = phi_alpha_old;
            phi_prime_a = phi_prime_alpha_old;
            b = alpha_new;
            phi_b = phi_alpha_new;
            [alpha_k,flag] = zoom(fval, gval,x_k,a,b,P_k,c1,c2,f0,gf0,phi_a,phi_b,phi_prime_a, varargin{:});
            return;
        end

        % Evaluate phi'(alpha_i)
        phi_prime_alpha_new = gval(x_k+alpha_new*P_k, varargin{:})'*P_k;

        if abs(phi_prime_alpha_new) <= -c2*phi_prime_0
            alpha_k = alpha_new;
            return;
        end

        if phi_prime_alpha_new >=0
            a = alpha_new;
            phi_a = phi_alpha_new;
            phi_prime_a = phi_prime_alpha_new;
            b = alpha_old;
            phi_b = phi_alpha_old;
            [alpha_k,flag] = zoom(fval, gval,x_k,a,b,P_k,c1,c2,f0,gf0,phi_a,phi_b,phi_prime_a, varargin{:});
            return;
        end

        dalpha = alpha_new - alpha_old;
        alpha_old = alpha_new;
        phi_alpha_old = phi_alpha_new;
        phi_prime_alpha_old = phi_prime_alpha_new;

        % choose alph_i \in [2alph_i - alpha_i_minus1,
        % min(alph_max,alph_i+tau1*(alpha_i-alph_i_minus1)]
        if alpha_max <= alpha_new+dalpha
            alpha_new = alpha_max;
        else
            alpha_new = (alpha_new+dalpha)/2 + min(alpha_max,alpha_new+tau1*dalpha)/2;
        end
    end
end %of wolfLineSearch

function [alpha_k,flag] =zoom(fval, gval,x_k,a,b,P_k,c1,c2,f0,gf0,phi_a,phi_b,phi_prime_a, varargin)
    % Sectioning
    flag = 0;
    epsilon = 1e-8;
    phi_prime_0 = gf0'*P_k;
    phi_0 = f0;
    counter = 0;
    % assume [[a_j,b_j]] is known and update for [[a_j+1,b_j+1]]
    while (1) %loop 2
        counter = counter + 1;
        % choosing appropriate alpha_j
        % using Quadratic interpolant of function phi(alpha)
        % using phi(a_j), phi(b_j), and Phi'(a_j) and then
        % minimizing of interpolant within interval [a_j,b_j]
        C1 = phi_a; C2 = (b-a)*phi_prime_a; C3 = phi_b -C1 - C2;
        q = @(z) C1 + C2*z + C3*z^3;
        if C3 ~= 0
            Z = -C2/(2*C3);
        else
            Z = -1; % a value outside of 0.1 and 0.5
        end
        if Z<0.1 || Z>0.5
            if q(0.1)<q(0.5)
                Z = 0.1;
            else
                Z = 0.5;
            end
        end

        alpha = a + Z*(b-a);

        % evaluate phi(alpha_j)
        phi_alpha = fval(x_k + alpha * P_k, varargin{:});

        % Break law is needed!!
        if (a-b)*phi_prime_a <= epsilon
            flag = 1; alpha_k = alpha;
            return;
        end

        if (phi_alpha > phi_0 + c1*alpha*phi_prime_0) || (phi_alpha >= phi_a)
            b = alpha;
            phi_b = phi_alpha;
        else
            % Evaluate phi'(alpha_j)
            phi_prime_alpha = gval(x_k + alpha*P_k, varargin{:})'*P_k;
            if abs(phi_prime_alpha) <= -c2*phi_prime_0
                alpha_k = alpha;
                return; % for loop 2
            end

            if (b-a)*phi_prime_alpha >= 0
                b = a;
                phi_b = phi_a;
            end

            a = alpha;
            phi_a = phi_alpha;
            phi_prime_a = phi_prime_alpha;

        end
        
        if counter > 100
            flag = 1; alpha_k = alpha;
            return;
        end
    end % of while
end % of zoom function