% Coding by Yang Hu


% Description: 
% To check the approximation ratio of 3.56, run the following code.
% To check the upper bound of the integrality gap 3.45, first comment lines 92, 136, 139, and then uncomment lines 91, 135, 138.

digits(32);
tic
function ret = solve_ratio(mu_l, mu_r, lprec, xprec, dl, dr, kl, kr, alpha_l, alpha_r)


    A = zeros(xprec + 1, lprec + 2);
    b = zeros(xprec + 1, 1);
    f = zeros(lprec + 2, 1);

    l = max(1, mu_l*alpha_l); % range of x: between mu_l*alpha_l and kr
    r = kr;

    delta = zeros(lprec, 1);
    lambda = zeros(lprec, 1);

    for i=1:lprec
        delta(i) = -(dl + (dr - dl) .* i ./ (lprec + 1));
        lambda(i) = log(1+delta(i));
    end

    for a=1:xprec+1
        for c=1:lprec
            x_l = l + (r - l) * (a - 1) / xprec;
            x_r = l + (r - l) * a / xprec;

            if a <= xprec
                % A(a,c) = -min(10000,exp(log(1 - delta(c)) * x_r + delta(c) * mu_l)); % trivial bound
                A(a,c) = -min(10000,exp(lambda(c) * x_r - ...
                    delta(c) * mu_l * (1-alpha_l) - (lambda(c)+0.0001) * alpha_l * mu_l)); % new bound
            else
                A(a,c) = 0;
            end
        end
    end

    for a=1:xprec+1
        x_l = l + (r - l) * (a - 1) / xprec;
        x_r = l + (r - l) * a / xprec;
        A(a,lprec + 1) = -1;
        if a <= xprec
            A(a,lprec + 2) = -max(0, kl-x_r);
        else
            A(a,lprec + 2) = 0;
        end
    end

    for c=1:lprec
        f(c) = 1;
    end
    f(lprec + 1) = 1;
    e = exp(1.0);
    f(lprec + 2) = ((e/(e-1))*(1-alpha_r)+(e/(e-1))^2*alpha_r) * mu_r;

    for a=1:xprec+1
        x_l = l + (r - l) * (a - 1) / xprec;
        b(a) = -max(0,log(kr)-log(x_l));
    end


    Aeq = [];
    beq = [];
    lb = zeros(lprec + 2, 1);
    ub = [];

    options = optimoptions(@linprog,'Algorithm','dual-simplex','Display','iter', 'TolCon', 1e-10);

    options.Display = 'off';

    [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub, options);
    ret = fval;
end

function ret = find(mu_l, mu_r, kl, kr, alpha_l, alpha_r, lprec, xprec)

    dl = 0;
    dr = 1;

    ret = solve_ratio(mu_l, mu_r, lprec, xprec, dl, dr, kl, kr, alpha_l, alpha_r);
end

function ret = search(mu_l, mu_r, kl, kr, alpha_l, alpha_r, step)
    ans1 = find(mu_l, mu_r, kl, kr, alpha_l, alpha_r, step * 2, step * 40);
    if step <= 15
        fprintf("search %d %d %d %d %d %d %d\n", mu_l, mu_r, kl, kr, alpha_l, alpha_r, ans1);
    end
    
    %limit = log(3.3); % check to see if any point in the range exceeds limit (for integrality gap)
    limit = log(3.55); % range exceeds limit for rounding algorithm
    if ans1 <= limit
        ret = ans1;
        return;
    end
    if mu_r - mu_l <= 0.0001 && kr - kl <= 0.0001 && alpha_r - alpha_l <= 0.0001
        ans1 = find(mu_l, mu_r, kl, kr, alpha_l, alpha_r, 300, 6000);
        if ans1 <= limit
            ret = ans1;
            return;
        end
        fprintf("dead %d %d %d %d %d %d %d\n", mu_l, mu_r, kl, kr, alpha_l, alpha_r, ans1);
        assert(false);
        return;
    else
        tans1 = 0;
        tans2 = 0;
        next_step = step + 1;
        [argvalue, argmax] = max([mu_r - mu_l, kr - kl, alpha_r - alpha_l]);
        if argmax == 1
            mid = (mu_l + mu_r) / 2;
            tans1 = search(mu_l, mid, kl, kr, alpha_l, alpha_r, next_step);
            tans2 = search(mid, mu_r, kl, kr, alpha_l, alpha_r, next_step);
        elseif argmax == 2
            mid = (kl + kr) / 2;
            tans1 = search(mu_l, mu_r, kl, mid, alpha_l, alpha_r, next_step);
            tans2 = search(mu_l, mu_r, mid, kr, alpha_l, alpha_r, next_step);
        else
            mid = (alpha_l + alpha_r) / 2;
            tans1 = search(mu_l, mu_r, kl, kr, alpha_l, mid, next_step);
            tans2 = search(mu_l, mu_r, kl, kr, mid, alpha_r, next_step);
        end
        ret = max(tans1, tans2);
        return;
    end
end

format longg
e = exp(1.0);



mu_l = 1;
% mu_r = 1000; % parameter of mu for the integrality gap
mu_r = 360000; % parameter of mu for the rounding algorithm
kl = 1;
% kr = 30000; % parameter of k for the integrality gap
kr = 20000000; % parameter of k for the rounding algorithm
alpha_l = 0;
alpha_r = 0;

realans = search(mu_l, mu_r, kl, kr, alpha_l, alpha_r, 1);

disp(realans);

toc
