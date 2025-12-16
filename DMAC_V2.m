function  [S, u_k, K, A_est, B_est] = DMAC_V2(S, PHI, y_k, r_k)
% y_k must be Nx1 vector  (x_{k+1} target for this update)
% u_k must be Mx1 vector  (u_k used in regressor)
% r_k must be sized of C*x
% S is the state of DMAC. Named S for the short name of State
% S MUST contain R0 and lambda

    % initialize DMAC
    if ~isfield(S, 'initFlag') || S.initFlag == 0
        S = initalizeDMAC(S, y_k, PHI);
        S.initFlag = 1;
    end

    Ny = S.length_y;  Nu = S.length_u;

    % update RLS
    [S.P_k, S.AB_est] = RLS_update(PHI, y_k, S.P_k, S.AB_est, S.lambda);
    % [S.AB_est] = gradDecent(S, PHI, y_k, S.AB_est);

    % extract A and B from AB_est
    A_est = S.AB_est(:, 1:Ny);  % where these are discrete estimations of A and B
    B_est = S.AB_est(:, (Ny+1):(Ny+Nu));

    [S, K] = getGainFSFI(S, A_est, B_est);
    [S, u_k] = getControlSignal(S,K,y_k,r_k);
    
    S.iteration = S.iteration + 1;
end

function [P_k_1, Theta_k_1, z_k] = RLS_update(PHI, y, P_k, Theta_k, lambda)
    % Recursive update for P_{k+1}
    PHI_trans_times_P_k = PHI' * P_k;

    gamma_k_inv = 1 / (lambda + PHI_trans_times_P_k * PHI);

    % covariance update
    P_k_1 = P_k/lambda - (1/lambda) * (P_k*PHI) * gamma_k_inv * (PHI_trans_times_P_k);

    % Recursive update for Theta_{k+1}
    Theta_k_1 = Theta_k + (y - Theta_k * PHI) * (PHI_trans_times_P_k * gamma_k_inv);

    z_k = norm(y - Theta_k * PHI);
end

function S = initalizeDMAC(S, y_k, PHI)
% checks if user defined FSFI and C, adds structure S.USE_FSFI_FLAG to S
    S.length_y = length(y_k);             % N (state dim)
    S.length_u = length(PHI(S.length_y+1:end));             % M (input dim)
    Ny = S.length_y;  Nu = S.length_u;

    S.P_k    = 1./S.R0(1) * eye(Ny + Nu);       % RLS covariance on [x;u]
    
    if ~isfield(S,'AB_est')
        S.AB_est = zeros(Ny, Ny + Nu);        % parameter matrix [A B]
    end

    [S.numberOfIntegralTerms, ~] = size(S.C);
    p = S.numberOfIntegralTerms;
    S.discreteErrorIntegral = zeros(p,1);
    
    % if S.UsePlaceAndLQR
        % require S.C (p×Ny)
        [pC, nC] = size(S.C);
        if pC ~= p || nC ~= Ny
            error('C does not have the correct dimensions (expected size of  y x x)') % number of integrator states x number of states
        end

        % LQR weights for augmented system
        if ~isfield(S,'Q')
            error('S.Q is not defined, weighting for LQR x states')
        elseif ~isfield(S,'R')
            error('S.R is not defined, weighting for LQR u states')
        end

        if ~isfield(S,'Z')
            % weighting for cross terms (S term in idare), size (Ny+p) x Nu
            S.Z = zeros(Ny + p, Nu);
        end
        if ~isfield(S,'E')
            % matrix, size (Ny+p) x (Ny+p)
            S.E = eye(Ny + p);
        end
    % else
    %     if ~isfield(S,'desiredPoles')
    %         error('S.desiredPoles is not defined')
    %     end
    % end
    if ~isfield(S,'K0')
        S.K = zeros(Nu, Ny + p); % gain sized for augmented state
    else
        S.K = S.K0;
    end

    % % grad decent initialization
    % S.aveGrad = [];
    % S.aveSqGrad = [];
    S.iteration = 1;
end

function [S, K] = getGainFSFI(S, A_est, B_est)
% Constructing Augmented State, assuming p = S.numberOfIntegralTerms integrals
    Ny = S.length_y;
    p  = S.numberOfIntegralTerms;
    
    [A_est_aug, B_est_aug] = augmentMatricies(A_est, B_est, S.C);

    % controllability check on augmented system
    if rank(ctrb(A_est_aug, B_est_aug)) ~= (Ny + p)
        % disp('Not Controllable!')
        K = S.K;
        S.acheivedPoles = nan(length(A_est_aug),1);
    else
        % if S.LQR
        %     %idare with Q (Ny+p)x(Ny+p), R (Nu x Nu), Z (Ny+p x Nu), E (Ny+p x Ny+p)
        %     % [~, K] = idare(A_est_aug, B_est_aug, S.Q, S.R, S.Z, S.E);
        %     % K = -K;
        %     % S.K = K;
        %     % S.acheivedPoles = eig(A_est_aug + B_est_aug*K);
        % else
        %     if S.iteration > 50
        %         K = -place(A_est_aug,B_est_aug, S.desiredPoles);
        %     else
        %         [~, K] = idare(A_est_aug, B_est_aug, S.Q, S.R, S.Z, S.E);
        %         K = -K;
        %     end
        % 
        %     S.acheivedPoles = eig(A_est_aug + B_est_aug*K); 
        %     S.K = K;     
        % end
        
        if S.UsePlaceAndLQR
            if S.iteration > 50
                    % K = -place(A_est_aug,B_est_aug, S.desiredPoles);
                    K = -placeAlex(A_est_aug, B_est_aug, S.desiredPoles);
            else
                    [~, K] = idare(A_est_aug, B_est_aug, S.Q, S.R, S.Z, S.E);
                    K = -K;
            end
        else
            [~, K] = idare(A_est_aug, B_est_aug, S.Q, S.R, S.Z, S.E);
            K = -K;
        end
            
        S.acheivedPoles = eig(A_est_aug + B_est_aug*K); 
        S.K = K;   
    end
end

function [S, u_k] = getControlSignal(S,K,y_k,r_k)
    % FSFI controller
    error_k = r_k - S.C*y_k;
    S.discreteErrorIntegral = S.discreteErrorIntegral + error_k;
    u_k = K*[y_k; S.discreteErrorIntegral];
end

function [K, acheivedPoles, pollErrors] = placeAlex(A,B,P)
% PLACEALEXFUN uses MIMO place algorithm from https://mech.novtex.ru/jour/article/download/22/368
% DOI: 10.17587/mau.19.11-18
% Use: u = K*x
% Closed Loop A matrix: A + B*K
    desiredPoles = P(:);
    n = size(A,1);
    % Helper: B_perp 
    computeB_perp = @(M) null(M','r')'; % (n-rank(M)) x n, orthonormal rows

    r0 = rank(B);
    L  = ceil(n / r0) - 1; % number of decomposition levels
    
    % Forward pass (Sec. 2, Eqs. 14–17)
    A_levels = cell(1,L+1); B_levels = cell(1,L+1); Bperp_levels = cell(1,L+1);
    A_levels{1} = A;  B_levels{1} = B;
    
    for k = 1:L
        A_prev = A_levels{k};
        B_prev = B_levels{k};
        B_perp = computeB_perp(B_prev); % B_perp * B_prev = 0, (paper's B_k-1^⊥)
        B_perp_pinv = pinv(B_perp); % Moore–Penrose
    
        % Level update:
        %   A_k = B_{k-1}^⊥ * A_{k-1} * B_{k-1}^{⊥+}
        %   B_k = B_{k-1}^⊥ * A_{k-1} * B_{k-1}
        A_k = B_perp * A_prev * B_perp_pinv;
        B_k = B_perp * A_prev * B_prev;
    
        A_levels{k+1} = A_k;
        B_levels{k+1} = B_k;
        Bperp_levels{k} = B_perp;                  % store B_{k-1}^⊥
    end
    
    A_L = A_levels{L+1};                           % final level matrices
    B_L = B_levels{L+1};
    
    % Build Phi_k blocks whose spectra equal desiredPoles
    Phi = cell(1,L+1);
    offset = 0;
    for k = 0:L-1
        Phi{k+1} = diag(desiredPoles(offset+1 : offset+r0));
        offset = offset + r0;
    end
    remDim = n - L*r0;
    Phi{L+1} = diag(desiredPoles(offset+1 : offset+remDim));
    
    % Backward pass (Sec. 3 & Fig. 1)
    F = cell(1,L+1);                    % F{1} = Ftheta, ..., F{L+1} = FL
   
    Bk_plus = pinv(B_L);
    F{L+1}  = Phi{L+1} * Bk_plus - Bk_plus * A_L;
    
    % For k = L-1 down to 0:
    for k = L-1:-1:0
        Ak = A_levels{k+1};
        Bk = B_levels{k+1};
        Bk_perp = Bperp_levels{k+1};    % this is B_k^⊥ (note the indexing shift)
    
        % B_k^- = B_k^+ − F_{k+1} * B_k^⊥
        Bk_minus = pinv(Bk) - F{k+2} * Bk_perp;
    
        % F_k = Phi_k * B_k^- − B_k^- * A_k
        F{k+1} = Phi{k+1} * Bk_minus - Bk_minus * Ak;
    end
    
    F0 = F{1}; % final state-feedback u = F0 * x
    
    % outputs
    K = F0; 
    acheivedPoles = eig(A + B * K); 
    pollErrors = acheivedPoles - desiredPoles;
end

% function [Theta_k] = gradDecent(S, PHI, y_k, Theta_k)
%     dJ_dtheta = -2 * (y_k - Theta_k * PHI) * PHI';
%     [Theta_k, S.aveGrad, S.aveSqGrad] = adamupdate_fast(Theta_k, dJ_dtheta, S.aveGrad, S.aveSqGrad, S.iteration, S.learnRate);
% end

function [Aa, Ba] = augmentMatricies(A, B, C)

    [ny, ~] = size(C);
    [~, nu] = size(B);
    nx = size(A,1);

    D = zeros(ny,nu);
 
    Aa = [A, zeros(nx,ny);
         -C,   eye(ny,ny)];
    Ba = [B; -D];
end