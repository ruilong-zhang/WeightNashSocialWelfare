digits(32);

function ret = solve_ratio(mu_l, mu_r, lprec, xprec, dl, dr)


    A = zeros(xprec + 1, lprec + 3);
    b = zeros(xprec + 1, 1);
    f = zeros(lprec + 3, 1);

    l = 1;
    r = mu_r * 5;

    delta = zeros(lprec, 1);

    for i=1:lprec
        delta(i) = dl + (dr - dl) .* i ./ (lprec + 1);
    end

    for a=1:xprec+1
        for c=1:lprec
            x_l = l + (r - l) * (a - 1) / xprec;
            x_r = l + (r - l) * a / xprec;

            if a <= xprec
                A(a,c) = -min(10000,exp(log(1 - delta(c)) * x_r + delta(c) * mu_l));
            else
                A(a,c) = 0;
            end
        end
    end

    for a=1:xprec+1
        x_l = l + (r - l) * (a - 1) / xprec;
        x_r = l + (r - l) * a / xprec;
        A(a,lprec + 1) = -1;
        A(a,lprec + 2) = -x_l;
        A(a,lprec + 3) = 1;
    end

    for c=1:lprec
        f(c) = 1;
    end
    f(lprec + 1) = 1;
    f(lprec + 2) = mu_r;
    f(lprec + 3) = -1;

    for a=1:xprec+1
        x_l = l + (r - l) * (a - 1) / xprec;
        b(a) = log(x_l);
    end


    Aeq = [];
    beq = [];
    lb = zeros(lprec + 3, 1);
    ub = [];

    options = optimoptions(@linprog,'Algorithm','dual-simplex','Display','iter', 'TolCon', 1e-10);

    [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub, options);
    ret = log(1+exp(1)/(exp(1)-1))+log(mu_r)+fval;
end

function ret = find(mu_l, mu_r)
    lprec = 30;
    xprec = 10000;

    dl = 0;
    dr = 1;

    ret = solve_ratio(mu_l, mu_r, lprec, xprec, dl, dr);
end

gap = 0.001;
x = 1:gap:4;
%gap = 0.1;
%x = 4:gap:300;
y = zeros(1,size(x,2));
v = 0;
for ii=1:size(x,2)
    y(1,ii) = find(x(ii),x(ii) + gap);
    v = max(v, y(1,ii));
end
figure
plot(x,y)