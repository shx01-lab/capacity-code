function [m_opt] = sollve_AB_SDP(M, theta, dA, dB)
% solve_sdp_max_m
%   maximize m
%   s.t. Tr(MF) >= 1 - theta
%        F >= 0
%        rho ⊗ I_B >= F
%        Tr_A(F) = I_B / m
%        Tr(rho) = 1
%        rho >= 0
%
% 通过替换 t = 1/m，把问题改写为：minimize t
% s.t. Tr(MF) >= 1 - theta
%      F >= 0
%      kron(rho, I_B) >= F
%      Tr_A(F) = t * I_B
%      Tr(rho) = 1, rho >= 0, t >= 0
%
% Inputs:
%   M     : (dA*dB) x (dA*dB) Hermitian matrix
%   theta : scalar
%   dA    : dim of subsystem A
%   dB    : dim of subsystem B
%
% Outputs:
%   m_opt : optimal m (approx), m_opt = 1/t_opt
%   t_opt : optimal t
%   rho_opt, F_opt : optimal variables
%   cvx_status_str : CVX status

nAB = dA * dB;
IB  = eye(dB);
M=(M+M')/2;
% 基本维度检查
assert(all(size(M) == [nAB, nAB]), 'M must be of size (dA*dB) x (dA*dB).');

cvx_begin sdp
    cvx_precision high

    variable rho(dA,dA) hermitian semidefinite
    variable F(nAB,nAB) hermitian semidefinite
    variable t nonnegative

    % --- Partial trace over A: Tr_A(F) ---
    % F is partitioned into dA x dA blocks, each block is dB x dB.
    expression TrA_F(dB,dB)
    TrA_F = 0;
    for i = 1:dA
        idx = (i-1)*dB + (1:dB);
        TrA_F = TrA_F + F(idx, idx);
    end

    % --- Constraints ---
    trace(rho) == 1;

    % Tr(MF) >= 1 - theta (取实部以避免数值导致的虚部极小误差)
    real(trace(M * F)) >= 1 - theta;

    % rho ⊗ I_B >= F  <=> kron(rho, I_B) - F is PSD
    kron(rho, IB) - F == semidefinite(nAB);

    % Tr_A(F) = t * I_B  (对应原式 I_B / m)
    TrA_F == t * IB;

    % --- Objective ---
    minimize(t)
cvx_end

t_opt = t;
cvx_status_str = cvx_status;

if t_opt <= 0
    m_opt = Inf;  % 理论上若 t=0 则 m 无穷大；通常不会发生
else
    m_opt = 1 / t_opt;
end

rho_opt = rho;
F_opt   = F;

end