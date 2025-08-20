classdef PURSE
    properties
        quadCon
        linCon
        centers
        radii
        kpts3d
        intrinsics
    end

    methods
        %% Constructor
        function obj = PURSE(mu,r,Y,P)
            % construct a purse from 2D balls
            % mu: 2 x N matrix, each column is a 2D point
            % r: 1 x N matrix, each column is the radius of a ball
            % Y: 3 x N matrix, each column is a 3D keypoint
            % P: 3 x 3 camera intrinsics, default P is identity
            if nargin < 4
                P = eye(3);
            end
            obj.intrinsics       = P;
            obj.centers          = mu;
            obj.radii            = r;
            obj.kpts3d           = Y;
            N                    = size(Y,2);
            A                    = zeros(12,12,N);
            b                    = zeros(12,N);
            for i = 1:N
                radius      = r(i);
                center      = mu(:,i);
                Lambda      = (1 / radius)^2 * eye(2);
                Yi          = Y(:,i);
                U           = [kron(Yi',P), P];
                u1          = U(1,:)';
                u2          = U(2,:)';
                u3          = U(3,:)';
                b(:,i)      = u3;
                Umu         = [u1-center(1)*u3, u2-center(2)*u3];
                A(:,:,i)    = Umu*Lambda*Umu' - u3*u3';
            end
            obj.quadCon = A;
            obj.linCon  = b;
        end

        %% normalize 2D keypoints (remove camera intrinsics)
        function y_norm = calibrate_kpts2d(obj,y)
            % y of shape (2,N)
            n_kpts = size(y,2);
            y      = [y;ones(1,n_kpts)];
            y_norm = obj.intrinsics \ y;
            y_norm = y_norm(1:2,:);
            % return y_norm also in shape (2,N)
        end
        
        %% Check membership
        function [isvalid,residuals] = checkMembership(obj,R,t)
            if(isrot(R) == false)
                isvalid = false;
                residuals = [];
            else
                N = size(obj.centers,2);
                s = [R(:);t];
                isvalid = true;
                all_residual1 = [];
                all_residual2 = [];
                for i = 1:N
                    bi = obj.linCon(:,i);
                    Ai = squeeze(obj.quadCon(:,:,i));
                    residual1 = [all_residual1,bi'*s];
                    residual2 = [all_residual2,s'*Ai*s];
                    if bi'*s <= 0
                        isvalid = false;
                        break
                    end
                    if s'*Ai*s > 0
                        isvalid = false;
                        break
                    end


                    residuals = [min(abs(residual1)),min(abs(residual2))];
                end
            end

        end

        %% RANSAG
        function [y,Y] = sample3pairs(obj, sampleMethod, edgeRatio)
            % sampleMethod:     "uniform"  : sample uniformly,
            %                   "weighted": sample according to weights, higher if the circle radius is smaller,
            %                   "minimum" : sample from min 3 candidates
            % edgeRatio: at what ratio will the sampled points be from the edge of the circle. 
            %            0 (default) means the sampled points will be uniformly sampled along the whole radius, i.e. [0, r]
            %            0.9 means the sampled points will be uniformly sampled along [0.9r, r]
            if ~exist("edgeRatio", "var")
                edgeRatio = 0;
            end
            assert(any(sampleMethod == ["uniform", "weighted", "minimum"]), "sampleMethod should be one of random, weighted, minimum, edge");
            if sampleMethod == "uniform"
                N   = size(obj.centers,2);
                ids = randsample(N, 3);
            elseif sampleMethod == "minimum"
                [~,ids] = mink(obj.radii,3);
                ids = ids';
            elseif sampleMethod == "weighted"
                ids = zeros(3, 1);
                weights = 1./obj.radii;
                for k = 1:3
                    new_rand = randsample(size(obj.centers, 2), 1, true, weights);
                    ids(k) = new_rand;
                    weights(new_rand) = 0;
                end
            end
            Y   = obj.kpts3d(:,ids);
            mu  = obj.centers(:,ids);
            r   = obj.radii(:,ids);
            v   = normc(randn(2,3));
            y   = mu + (r .* (1-rand(1,3)*(1-edgeRatio))) .* v;
        end

        function [Rset,tset,Ravg,tavg] = ransag(obj,T,sampleMethod, edgeRatio)
            % T: maximum number of iterations
            % sampleMethod:     "uniform"  : sample uniformly,
            %                   "weighted": sample according to weights, higher if the circle radius is smaller,
            %                   "minimum" : sample from min 3 candidates
            % edgeRatio: at what ratio will the sampled points be from the edge of the circle. 
            %            0 (default) means the sampled points will be uniformly sampled along the whole radius, i.e. [0, r]
            %            0.9 means the sampled points will be uniformly sampled along [0.9r, r]
            if ~exist("edgeRatio", "var")
                edgeRatio = 0;
            end
            Rset = {};
            tset = {};
            Ravg = eye(3);
            tavg = zeros(3,1);
            Rsum = zeros(3,3);
            tsum = zeros(3,1);

            % flag_fmincon = 0;
            % 
            % if(sampleMethod == "edge_fmincon")
            %     sampleMethod = "uniform";
            %     flag_fmincon = 1;
            % end
            for i = 1:T
                [y,Y] = obj.sample3pairs(sampleMethod, edgeRatio);
                % remove camera intrinsics from 2D keypoints
                y     = obj.calibrate_kpts2d(y);
                [R,t] = solveP3P(y,Y);
                if ~isempty(R)
                    for j = 1:size(t,2)
                        if(isrot(R(:,:,j)) == false)
                            continue;
                        end
                        inside = obj.checkMembership(R(:,:,j),t(:,j));
                        if inside
                            Rsum = Rsum + R(:,:,j);
                            tsum = tsum + t(:,j);
                            Rset = [Rset;{R(:,:,j)}];
                            tset = [tset;{t(:,j)}];
                        end
                    end
                end
            end

            n = length(Rset);
            
            if any(isnan(Rsum)) | any(isinf(Rsum))
                idx1 = find(isnan(Rsum));
                idx2 = find(isinf(Rsum));

                Rsum(idx1) = 0;
                Rsum(idx2) = 0;
            end
            
            if n > 0
                tavg = tsum / n;
                Ravg = project2SO3(Rsum);
            end
        end

        function [R_center, R_trajectory, R_targets, t_targets, R_upperbounds, R_lowerbounds, R_gaps] = iterateRotationCenter(obj, max_iter, R_center_init)

            if ~exist("R_center_init", "var")
                [R_set, t_set, R_avg, t_avg] = obj.ransag(2000, "uniform");
                R_center = R_set{1};
            else
                R_center = R_center_init;
            end
            
            R_trajectory = zeros(3,3,max_iter);
            R_targets = zeros(3,3,max_iter);
            t_targets = zeros(3,max_iter);
            R_upperbounds = zeros(max_iter, 1);
            R_lowerbounds = zeros(max_iter, 1);
            R_gaps = zeros(max_iter, 1);
            start_time = tic;

            for t = 1:max_iter
                fprintf("\n----------------- Start iteration %d -----------------\n", t);
                iter_start_time = tic;
                step = 1/t;
                [temp_R, temp_t, info] = maxRotationDist(R_center, obj, 2);
                R_targets(:, :, t) = temp_R;
                t_targets(:, t) = temp_t;
                R_center = Geodesic(R_center, temp_R, step);
                R_trajectory(:,:,t) = R_center;
                R_upperbounds(t) = rad2deg(acos(1 - info.f_sdp/4));
                R_lowerbounds(t) = rad2deg(acos(1 - info.f_est/4));
                R_gaps(t) = rad2deg(acos(1 - info.gap/4));
                fprintf("time: %.3f\n", toc(iter_start_time));
            end
            fprintf("Total time: %.3f, avg time: %.3f\n", toc(start_time), toc(start_time)/max_iter);
        end

        function normalize(obj)
            N = size(obj.centers,2);
            quadCon = zeros(12,12,N);
            linCon = zeros(12,N);
            for i = 1:N
                Ai = squeeze(obj.quadCon(:,:,i));
                Bi = squeeze(obj.linCon(:,i));
                quadCon(:,:,i) = Ai / norm(Ai(:));
                linCon(:,i) = Bi / norm(Bi(:));
            end
            obj.quadCon = quadCon;
            obj.linCon = linCon;
        end

    end


end

function c_out = Geodesic(input, target, ratio)

    cin = rotm2quat(input);
    qf = rotm2quat(target);
    qua_cin = quaternion(cin);
    qua_f = quaternion(qf);
    qua_cout = slerp(qua_cin, qua_f, ratio);
    c_out = quat2rotm(qua_cout);
end