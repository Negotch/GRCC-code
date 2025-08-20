function sol = GRCC(problem, kappa, P, Q, d)
    x = problem.vars;
    problem.objective = x'*x; %Place holder, doesn't matter.
    [SDP,~] = dense_sdp_relax(problem,kappa);
    sol = general_lass_sedumi(SDP,P,Q,d);
end