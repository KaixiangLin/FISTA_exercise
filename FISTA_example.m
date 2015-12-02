% optimize the objective: min_x 0.5 || Ax - y ||^2 + rho1 * ||x||_1
% A: n by d. 
% y: n by 1.

function [x, funcVal] = FISTA_example(A, y, rho1, maxIter, tol)

dimension = size(A, 2);
funcVal = [];

% flags of optimization algorithms . 
bFlag=0; % this flag tests whether the gradient step only changes a little
tFlag=3; % the termination criteria

% initialize a starting point
x0 = randn(dimension, 1);
xk   = x0; % x_k     
xk_1 = x0; % x_{k-1}

tk   = 1; % t_k
tk_1 = 0; % t_{k-1}

% initialize other variables. 
L    = 0.1; % initial L       << here we use 0.1 so you can see how line search goes. 
eta  = 1.1; % increment of L 

iter = 0;
Lk = L;
while iter < maxIter
    alpha = (tk_1 - 1) /tk;
    
    % get current search point by a linear combination 
    yk = (1 + alpha) * xk - alpha * xk_1;
    
    % compute function value and gradients of the search point
    fg_yk  = gradVal_eval (yk); % gradient 
    fv_yk  = funVal_eval  (yk); % function value
    
    
    % start line search 
    while true
        pL_yk = l1_projection(yk-1/Lk*fg_yk,rho1/Lk);
        Fv_plyk = funVal_eval(pL_yk)+ rho1 * sum(abs(pL_yk));
        q_apro = quadratic_approx(pL_yk,yk,Lk);
        % the line search procedure goes here!
        if (Fv_plyk   <= q_apro)
            break;
        end
        Lk = Lk * eta;
        
    end
    
 
    % update current and previous solution.
    xk_1 = xk;
    xk = pL_yk;
    
    % concatenate function value (smooth part + non-smooth part).
 
    funcVal = cat(1, funcVal, Fv_plyk);
    
    if (bFlag)
        % fprintf('\n The program terminates as the gradient step changes the solution very small.');
        break;
    end
    
    % test stop condition.
    switch(tFlag)
        case 0
            if iter>=2
                if (abs( funcVal(end) - funcVal(end-1) ) <= tol)
                    break;
                end
            end
        case 1
            if iter>=2
                if (abs( funcVal(end) - funcVal(end-1) ) <= tol* funcVal(end-1))
                    break;
                end
            end
        case 2
            if ( funcVal(end)<= tol)
                break;
            end
        case 3
            if iter>= maxIter
                break;
            end
    end
    
    % update other variables. 
    iter = iter + 1;
    tk_1 = tk;
    tk = 0.5 * (1 + (1+ 4 * tk^2)^0.5);
    
end

x = pL_yk;


% private functions

    
    % gradient of smooth part f at a given point x. 
    function [grad_x] = gradVal_eval(x)
        % gradient goes here
%         error('TO BE IMPLEMENTED')
         grad_x = A'*(A*x - y);
    end

    % function value of smooth part f.
    function [funcVal] = funVal_eval (x) 
        % function value goes here. 
        funcVal = 0.5*sum((A*x - y).^2);
    end

    % Quadratic approximation
    function q = quadratic_approx(pl_yk, yk, L)
        % this function calculate the quadratic approximation of function. 
        q = funVal_eval(yk) + (pl_yk - yk)'* gradVal_eval(yk) + L/2*sum((pl_yk-yk).^2) + rho1 * sum(abs(pL_yk));
    end 

end

% projection 
function z = l1_projection (v, beta)
    % this projection calculates
    % argmin_z = 0.5 * \|z-v\|_2^2 + beta \|z\|_1
    % z: solution
    % l1_comp_val: value of l1 component (\|z\|_1)
    z = sign(v) .* max(0, abs(v)- beta); 
end

