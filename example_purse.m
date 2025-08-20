
close all, clear all, clc

load("purse_data.mat")
addpath("./utils")
addpath("p3p_problem");
addpath("spotless/")


%% GRCC for Translation with kappa = 2
kappa = 2;

P = [zeros(3,9),eye(3)];
Q = eye(3);

sol = modified_purse_lass_sedumi(purse,P,Q,kappa);

translation_rad_cheb = sol.upper_bound;

%% GRCC for Rotation with kappa = 3

kappa = 3;

[R_set,t_set,Ravg,tavg] = purse.ransag(50000, "uniform");
R_set = cat(3, R_set{:});
% check the max distance
[~,~,avg_info] = maxRotationDist(Ravg,purse,2);
rotation_rad_avg = rad2deg(acos(1 - avg_info.f_sdp / 4));
if rotation_rad_avg > 90
    fprintf("The distance is too large: %.2f, skip this one!\n",rad_avg)
end
[~,~,t_info] = maxTranslationDist(tavg,purse,2);
translation_rad_avg = sqrt(t_info.f_sdp);

%% Calculate the performance of the average pose
q_avg = rotm2quat(Ravg);
q_avg = q_avg(:);

%% Calculate the performace of the 'best‘
[center,sol] = get_quat_cent2(purse,kappa,q_avg);
center = center/norm(center);
r_cent = quat2rotm(center');
[R2,t2,Rinfo2] = maxRotationDist(r_cent,purse,2);
rotation_rad_cheb = rad2deg(acos(1 - Rinfo2.f_sdp / 4));

fprintf("\n Translation:\n Avg rad: %.2f\n Cheb rad: %.2f \n", translation_rad_avg, translation_rad_cheb);
fprintf("\n Rotation:\n Avg rad: %.2f\n Cheb rad: %.2f \n", rotation_rad_avg, rotation_rad_cheb);


function [center,sol] = get_quat_cent2(purse,kappa,q_avg)
    [g,h,x] = get_ineq(purse);
    g = [g;-q_avg'*x(1:4)];
    f = -x'*x;    problem.vars = x;
    % Need to renormalize 
    for i = 1: length(g)
         gi = g(i);
         g(i) = subs(gi, x(7), 20*x(7));
     end
     for j = 1: length(h)
         hj = h(j);
         h(j) = subs(hj, x(7), 20*x(7));
     end
    [g, ~] = msspoly_clean(g, x, 1e-12, true);
    [h, ~] = msspoly_clean(h, x, 1e-12, true);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    problem.objective = f;
    problem.inequality = g;
    problem.equality = h;
    [SDP,~] = dense_sdp_relax(problem,kappa);
    d = 7;
    P = [eye(4),zeros(4,3)];
    Q = eye(4);
    sol = general_lass_sedumi(SDP,P,Q,d);
    center = sol.a;
end
    
    
function [g,h,vars,info] = get_ineq(purse)
    A = purse.quadCon;
    B = purse.linCon;
    d       = 7; 
    x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
    r       = x(1:4);
    t       = x(5:7);
    R       = [2*(r(1)^2+r(2)^2)-1, 2*(r(2)*r(3)-r(1)*r(4)), 2*(r(2)*r(4) + r(1)*r(3));...
                2*(r(2)*r(3)+r(1)*r(4)), 2*(r(1)^2+r(3)^2)-1, 2*(r(3)*r(4) - r(1)*r(2));...
                2*(r(2)*r(4)-r(1)*r(3)), 2*(r(3)*r(4)+r(1)*r(2)), 2*(r(1)^2+r(4)^2)-1;];
    h       = r'*r - 1;
    y = [R(:);t]; 
    g       = [30^2 - t'*t];
    for i = 1:size(A,3)
        Ai = squeeze(A(:,:,i));
        g  = [g; y' * (-Ai) * y];
    end
    slack = 0;
    for i = 1:size(B,2)
        Bi = squeeze(B(:,i));
        Bi = Bi(:);
        Bi = Bi/norm(Bi,'fro');
        g  = [g; y' * Bi - slack];
    end
    vars = x;
    info.R = R;
end

function [cleaned_pop, scale_vec] = msspoly_clean(pop, v, tol, if_scale)
    % pop: msspoly polynomial vector of size (N, 1)
    % v: msspoly variable of size (n, 1)
    % tol: numerical tolerance
    % if_scale: scale each coefficient to [-1, 1] or not
    % scale_vec: contain the max absolute value of coefficients in each polynomial, of size (N, 1)
    N = size(pop, 1);
    pop = [pop; v'*v];
    [~, degmat, coefmat] = decomp(pop); % degmat of size (N_term, n), coefmat of size (N, N_term)

    cleaned_pop = [];
    scale_vec = [];
    for k = 1: N
        p = pop(k);
        p_coef = coefmat(k, :);
        [row, ~, val] = find(p_coef');
        max_abs_val = max(abs(val));                
        if if_scale
            val = val / max_abs_val;
            scale_vec = [scale_vec; max_abs_val];
        else
            scale_vec = [scale_vec; 1];
        end
        cleaned_p = 0;
        for i = 1: length(row)
            term_degmat = degmat(row(i), :);
            [row_deg, ~, val_deg] = find(term_degmat');
            term = 1;
            for j = 1: length(row_deg)
                var_id = row_deg(j);
                term = term * v(var_id)^val_deg(j);
            end
            if abs(val(i)) > tol
                cleaned_p = cleaned_p + val(i) * term;
            end
        end
        cleaned_pop = [cleaned_pop; cleaned_p];
    end
end