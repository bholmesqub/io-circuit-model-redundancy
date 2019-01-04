function [d, D, alpha ] = detectRedundancy(theta, kappa, varargin)
% detectRedundancy  A function to determine if there are redundant
% parameters in a model given a parameter vector and exhaustive summary

% Author:   Ben Holmes
% Date:     2019/01/04

% Inputs:
%   - theta:    Parameter vector (symbolic)
%   - kappa:    Exhaustive summary of model (symbolic).
%   - simplify: 'log', 'log-diag', or 'none'. Can drastically improve
%   computation time of 

% Outputs:
%   - d:        Number of redundant parameters
%   - D:        Derivative matrix of dkappa / dtheta (symbolic).
%   - alpha:    Nullspace of derivative matrix D (symbolic). Output is left
%   as optional as it can be incredibly computationally expensive to
%   calculate the rank.

% References
% "Determining the parametric structure of models" - D. J. Cole et. al.,
% 2010, Mathematical Biosciences
% All of Diana Cole's Maple code

%% Input parsing

p = inputParser;

% Theta and kappa must be symbolic vectors
isSymVec =@(x) strcmp(class(x),'sym') && any(size(x) == 1);

p.addRequired('theta',isSymVec);
p.addRequired('kappa',isSymVec);

% Simplify is a character array that must be log log-diag or none
p.addOptional('simplify','log-diag',@(x) isa(x,'char')...
                                    &&(strcmp(x,'log')...
                                    || strcmp(x,'log-diag')...
                                    || strcmp(x,'none')));

% Parse
p.parse(theta, kappa, varargin{:});
simplify = p.Results.simplify;

%% Find jacobian matrix for redundancy analysis

% If required find the full jacobian matrix
if nargout >= 2
    D = jacobian(kappa, theta);
end

% If simplification is indicated
% Apply log function to exhaustive summary
if strcmp(simplify(1:3),'log')
    Dsimp = jacobian(log(kappa), theta);
    
    % Apply further simplification by multiplication by diagonal matrix of
    % theta.
    if strcmp(simplify, 'log-diag')
        Dsimp = (diag(theta)*(Dsimp.')).';
    end
elseif strcmp(simplify, 'none')
    if ~exist('D','var')
        Dsimp = jacobian(kappa, theta);
    else
        Dsimp = D;
    end
end

%% Calculate the number of redundant parameters.

nParams = length(theta);
nEstimable = rank(Dsimp);

d = nParams - nEstimable;

%% Report and calculate the alpha vector for finding estimable combinations

alpha = [];
if d == 0
    fprintf('\n\tNo redundant parameters: all model parameters are locally estimable.\n\n');
elseif d>0
    fprintf('\n\tRedundancy detected: %d redundant parameter(s), %d locally estimable combination(s)\n\n', d, nParams-d);
    if nargout >= 3
        alpha = null(D);
    end
end

end