
% GKr = 0.035 * sqrt(Ko / 5.4);

% to test: 
% state= [1,0,0,0,0,0,0]
 % para.dGkr_PKA=0
 % para.Vkr_PKA=0
% [ydot, IKr_out] = IKr_ODE(state, para, 0, -1, 0.1, -86, 5.4)

function [ydot, IKr_out] = IKr_ODE(state, Vkr_PKA, dGkr_PKA, Drug_D, Vm, GKr, Ek, Ko)
    C1     = state(1);
    C2     = state(2);
    C3     = state(3);
    I      = state(4);
    O      = state(5);
    I_star = state(6);
    O_star = state(7);

%     alpha = 1.8e-1 * exp(4.6e-2 * (Vm - 2.6e1 + para.Vkr_PKA));
%     beta = 1.9e-1 * exp(-4.6e-2 * (Vm + para.Vkr_PKA));
        alpha = 1.8e-1 * exp(4.6e-2 * (Vm - 2.6e1 + Vkr_PKA));
    beta = 1.9e-1 * exp(-4.6e-2 * (Vm + Vkr_PKA));
    alpha1 = 2.7e0;
    beta1 = 1.6e0;
%     alpha2 = 7.4e-2 * exp(1.8e-2 * (Vm - 7.6e1 + para.Vkr_PKA));
%     beta2 = 2.7e-3 * exp(-6.5e-3 * (Vm + para.Vkr_PKA));
        alpha2 = 7.4e-2 * exp(1.8e-2 * (Vm - 7.6e1 + Vkr_PKA));
    beta2 = 2.7e-3 * exp(-6.5e-3 * (Vm + Vkr_PKA));
    alphai = 1.7e-1 * exp(-1.7e-2 * (Vm + 1.6e1) * (4.5 / Ko));
    betai = 8.5e-1 * exp(1.8e-2 * Vm * (4.5 / Ko)^0.3);
    mu = alphai * beta2 / betai;

    Ka = 2.01e1 / 1e3; % (1/uM/ms)
    lA = 2.87e-4 / 1e3; % (1/uM/ms)
    Ki = 7.14e-1 / 1e3; % (1/uM/ms)
    lI = 8.99e-6 / 1e3; % (1/uM/ms)

    dC1     = beta * C2 - alpha * C1;
    dC2     = alpha * C1 + beta1 * C3 - (beta + alpha1) * C2;
    dC3     = alpha1 * C2 + beta2 * O + mu * I - (beta1 + alpha2 + alpha2) * C3;
    dI      = alpha2 * C3 + betai * O + lI * I_star - (mu + alphai + Ki * Drug_D) * I;
    dO      = alpha2 * C3 + alphai * I + lA * O_star - (beta2 + betai + Ka * Drug_D) * O;
    dI_star = Ki * Drug_D * I - lI * I_star;
    dO_star = Ka * Drug_D * O - lA * O_star;

    ydot = [dC1; dC2; dC3; dI; dO; dI_star; dO_star];

    % Calculate IKr
%     IKr_out = (1 + para.dGkr_PKA) * GKr * O * (Vm - Ek);
        IKr_out = (1 + dGkr_PKA) * GKr * O * (Vm - Ek);

end
