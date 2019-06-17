%%% Find parameter redundancy in the transfer function of the RC circuit.

clear;

%% Parameters

% Parameters
syms R1 C1

% Parameter vector
theta = [R1 C1];

%% Find transfer function

N_R = [1 -1]; % Resistor connected between node 1 and 2
N_C = [0  1]; % Capacitor connected between node 2 and ground
N_V = [1  0]; % Voltage source connected between node 1 and ground
N_O = [0  1]; % Output voltage measured between node 2 and ground

syms s  % Laplace variable

% Conductances
G_R = 1/R1;
G_C = s*C1;

% System matrix
S = [N_R.'*G_R*N_R + N_C.'*G_C*N_C N_V';...
        N_V, 0];
    
% Vector containing no current sources and 1 voltage source
U = [0; 0; 1];

% All transfer functions (Vi->Vi, Vi->V2, Vi->Iv)
rcTransfers = S\U;

% H(s) = 1/(1 + s*C1*R1)
ioTransfer = [N_O 0]*rcTransfers;

%% Extract exhaustive summary

% Find numerator and denominator of transfer function
[num, den] = numden(ioTransfer);

% Find all coefficients of transfer function
kappaInit = [coeffs(num, s), coeffs(den, s)];

% Remove irrelevant elements
kappa = unique(kappaInit(has(kappaInit, theta)));

%% Test transfer function

% Find the first order derivatives of kappa vs theta
J_kappa = jacobian(kappa, theta);

nParameters = length(theta);

nEstimable = rank(J_kappa);

% Report number of redundant parameters.
nRedundant = nParameters - nEstimable;
fprintf('%d redundant parameters.', nRedundant);