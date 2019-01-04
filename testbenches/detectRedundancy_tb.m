%%% Test bench to validate the behaviour of detectRedundancy
% Date: 2019/01/04
% Author: Ben Holmes

clear;
clc;

syms Gr1 Gr2 Gr3

% Parameters
theta = [Gr1 Gr2 Gr3];

% Combinations that could be found in a model
kappaRedundant = [Gr1*Gr2 Gr1*Gr2*Gr3 Gr1*Gr2*Gr3^2];
kappaEstimable = [Gr1*Gr2 Gr1*Gr3 Gr2*Gr3];

% Simplification methods to improve computation time
simplificationMethods = {'none','log','log-diag'};
nMethods = length(simplificationMethods);

%% Test detectRedundancy for each simplification

for nn=1:nMethods
    %% Test that redundancy is detected when it is present

    [d, D, alpha] = detectRedundancy(theta, kappaRedundant, simplificationMethods{nn});

    if d ~= 1
        error('kappaRedundant only has two estimable combinations.');
    end
    
    if any(any(jacobian(kappaRedundant, theta) ~= D))
        error('D must be jacobian of kappa relative to theta.');
    end
    
    if isempty(alpha)
        disp(null(D));
        error('Alpha should contain nullspace of D.');
    end

    %% Test that estimable combinations are correctly detected

    [d, D, alpha] = detectRedundancy(theta, kappaEstimable, simplificationMethods{nn});

    if d ~= 0
        error('kappEstimable should have 3 estimable parameters.');
    end
    
    if any(any(jacobian(kappaEstimable, theta) ~= D))
        error('D must be jacobian of kappa relative to theta.');
    end
    
    if ~isempty(alpha)
        disp(null(D));
        error('Nullspace of D should be undefined.');
    end
    
end

%% Test for nonsense inputs

% Check for non-symbolic theta
try
    d = detectRedundancy([1 2 3], kappaEstimable, 'none');
catch exception
    disp(exception.message);
end

% Check for non-symbolic kappa
try
    d = detectRedundancy(theta, [1 2 3], 'none');
catch exception
    disp(exception.message);
end

% Check for matrix theta
try
    d = detectRedundancy([Gr1 Gr2; Gr3 Gr3], kappaEstimable, 'none');
catch exception
    disp(exception.message);
end

% Check for matrix kappa
try
    d = detectRedundancy(theta, [Gr1 Gr2; Gr3 Gr3], 'none');
catch exception
    disp(exception.message);
end

% Check for non-string simplify
try
    d = detectRedundancy(theta, kappaEstimable, 1);
catch exception
    disp(exception.message);
end

% Check for wrong string simplify
try
    d = detectRedundancy(theta, kappaEstimable, 'wrong');
catch exception
    disp(exception.message);
end

% Check for too many inputs
try
    d = detectRedundancy(theta, kappaEstimable, 'log', 1);
catch exception
    disp(exception.message);
end

% Check for too few inputs
try
    d = detectRedundancy(theta);
catch exception
    disp(exception.message);
end