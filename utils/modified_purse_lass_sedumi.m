function sol = modified_purse_lass_sedumi(purse,P,Q,kappa)
    %% Generate relax problem:
    A = purse.quadCon;
    B = purse.linCon;
    d       = 12; 
    x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
    r       = x(1:9);
    t       = x(10:12);
    R       = reshape(r,3,3);
    c1      = R(:,1);
    c2      = R(:,2);
    c3      = R(:,3);
    f       = - x'*x;
    h       = [c1'*c1 - 1;
               c2'*c2 - 1;
               c3'*c3 - 1;
               c1'*c2;
               c2'*c3;
               c3'*c1;
               cross(c1,c2) - c3;
               cross(c2,c3) - c1;
               cross(c3,c1) - c2];
           
    % g       = [10^2 - sum(t.^2)];
    % g = [];
    g       = [10^4 - sum(t.^2)];
    for i = 1:size(A,3)
    Ai = squeeze(A(:,:,i));
    g  = [g; x' * (-Ai) * x];
    end
    slack = 1e-3;
    for i = 1:size(B,2)
        Bi = squeeze(B(:,i));
        Bi = Bi(:);
        g  = [g; x' * Bi - slack];
    end
    problem.vars = x;
    problem.objective = f;
    problem.inequality = g;
    problem.equality = h;

    sol = GRCC(problem, kappa, P, Q, d);
end