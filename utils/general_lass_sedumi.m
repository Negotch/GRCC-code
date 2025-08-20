function sol = general_lass_sedumi(SDP,S,Q,d)
    nm = size(S,1);
    assert(size(S,2) == d);
    assert(size(Q,1) == nm);
    assert(size(Q,2) == nm);
    Q_inv = inv(Q);
    T = S'*Q*S;
    %% modify blk, At, b, and C
    blk = SDP.blk;
    At  = SDP.At;
    b   = SDP.b;
    C   = SDP.C;
    sedumi_data = SDP.sedumi;
    s = SDP.sedumi.K.s;
    s_At = SDP.sedumi.At;

    % Add Semidefinite variable S to blk
    % [I, a
    %  a', t] \succeq 0
    blk = [blk;{'s',nm+1}];
    s = [s;nm+1];
    sedumi_data.K.s = s;
    n_blks = size(blk,1);

    % first, append all zeros to S's first m constraints
    At = [At;
          {sparse(triangle_number(blk{end,2}),length(b))}];

    s_At = [s_At;sparse(s(end)^2,length(b))];

    % second, add constraint the top-left dxd block of S is identity
    S_eye_con_lhs = {};
    S_eye_con_rhs = [];
    for i = 1:nm
        for j = i:nm
            Aij = sparse([i],[j],[1],blk{end,2},blk{end,2});
            bij = Q_inv(i,j);
            if i ~= j
                Aij = (Aij + Aij')/2; % make it symmetric if not
                % bij = 0; % off-diagonal is zero
            end
            S_eye_con_lhs = [S_eye_con_lhs;{Aij}];
            S_eye_con_rhs = [S_eye_con_rhs;bij];
        end
    end

    At{end} = [At{end}, sparsesvec(blk(end,:), S_eye_con_lhs)];
    b       = [b;sparse(S_eye_con_rhs)];
    for k = 1:n_blks-1
        At{k} = [At{k}, sparse(triangle_number(blk{k,2}), length(S_eye_con_rhs))];
    end

    s_At_add_2 = [sparse(size(s_At,1) - s(end)^2,length(S_eye_con_rhs));sparsevec(blk(end,:),S_eye_con_lhs)];
    s_At = [s_At,s_At_add_2];

    % third, add constraint the right blk is equal to a
    S_a_con_lhs = {};
    mom_a_con_lhs = {};
    S_a_con_rhs = sparse(nm,1); % all zeros
    for i = 1:nm
        Ai_S = sparse([i],[nm+1],[-1],blk{end,2},blk{end,2}); % constraint to be added to S
        % Ai_mom = sparse([i+1],[1],[-1],blk{1,2},blk{1,2}); % constraint to be added to moment matrix
        Ai_mom = sparse(2:d+1,[1],S(i,:),blk{1,2},blk{1,2}); % constraint to be added to moment matrix
        % make symmetric
        Ai_S = (Ai_S + Ai_S');
        Ai_mom = (Ai_mom + Ai_mom');

        S_a_con_lhs = [S_a_con_lhs; {Ai_S}];
        mom_a_con_lhs = [mom_a_con_lhs; {Ai_mom}];
    end

    At{1} = [At{1}, sparsesvec(blk(1,:), mom_a_con_lhs)];
    At{end} = [At{end}, sparsesvec(blk(end,:), S_a_con_lhs)];
    for k = 2:n_blks-1
        At{k} = [At{k}, sparse(triangle_number(blk{k,2}),nm)];
    end

    b = [b;S_a_con_rhs];

    s_At_add_3 = [sparsevec(blk(1,:), mom_a_con_lhs);sparse(size(s_At,1) - s(end)^2 - (blk{1,2})^2,nm);sparsevec(blk(end,:), S_a_con_lhs)];
    s_At = [s_At,s_At_add_3];

    % lastly, modify the cost function
    % the cost function is t - (a_1^2 + ... + a_d^2)
    idxs_1 = repmat(1:d,1,d);
    idxs_2 = [];
    for i = 1:d
        idxs_2 = [idxs_2,i*ones(1,d)];
    end
    % C1 = sparse(2:d+1,2:d+1,-1*ones(1,d),blk{1,2},blk{1,2});
    C1 = sparse(idxs_1+1,idxs_2+1,-T(:),blk{1,2},blk{1,2});
    Cend = sparse([nm+1],[nm+1],[1],blk{end,2},blk{end,2});
    C{1} = C1;
    C    = [C;{Cend}];

    sc          = [];
    for i = 1:length(C)
        sc      = [sc;sparsevec(blk(i,:),C(i))];
    end
    s_c    = sc;


    % re-assigne
    SDP.blk = blk;
    SDP.At  = At;
    SDP.b   = b;
    SDP.C   = C;

    SDP.sedumi.K.s = s;
    SDP.sedumi.At = s_At;
    SDP.sedumi.c = s_c;
    SDP.sedumi.b = b;

    % convert to sedumi
    % [At,new_b,c,K] = SDPT3data_SEDUMIdata(SDP.blk,SDP.At,SDP.C,SDP.b);
    % sedumi.At = At;
    % sedumi.b  = b;
    % sedumi.c  = c;
    % sedumi.K  = K;


    % convert to MOSEK
    % prob = convert_sedumi2mosek(At, new_b, c, K);
    prob = convert_sedumi2mosek(s_At, b, s_c, SDP.sedumi.K);


    [~,res]    = mosekopt('minimize info',prob);
    % [~,res]    = mosekopt('minimize log(mosek_log.txt)',prob);
    % [~,res]    = mosekopt('minimize echo(0)',prob);
    
    
    [Xopt,~,~,obj] = recover_mosek_sol_blk(res,SDP.blk);
    
    sol.a = S*Xopt{1}(2:d+1,1);
    sol.upper_bound = sqrt(-obj(1));
    sol.Xopt = Xopt{1};
    % sol.flag = purse.checkMembership(R,sol.a);
    % disp(sol.flag)
end