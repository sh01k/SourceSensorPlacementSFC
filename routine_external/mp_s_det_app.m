function [zhat, L, ztilde, Utilde] = mp_s_det_app(G, K)
% Solves the problem
%     maximize log det (sum_{m=1}^M z_m g_m g_m^T) + kappa sum_{m=1}^m(log(z_m)+ log(1-z_m))
%	  subject to sum(z) = K
%        	     0 <= z_m <= 1, m=1,..., M
%     http://ieeexplore.ieee.org/document/4663892/
% INPUT
%     G                Measurement matrix
%     K                Number of sensors to be selected
% OUTPUT
%     zhat             Estimate z of binary
%     L                Objective value
%     ztilde           Estimate z of real number 
%     Utilde           Upper bound of optimization problem
% 
% Modified from original code available here: 
%     www.stanford.edu/~boyd/papers/sensor_selection.html
% 
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

% Newton's method parameters
MAXITER  = 30;
NT_TOL = 1e-3;
GAP = 1.005;
% Backtracking line search parameters
alpha = 0.01;
beta = 0.5;

[M, N] = size(G);
z = ones(M,1)*(K/M); % initialize
% g = zeros(M,1);
ones_M = ones(M,1);
kappa = log(GAP)*N/M; 
% guarantees GM of lengths of semi-axes of ellipsoid corresponding to 
% ztilde <= 1.01 from optimal

fprintf('\nIter.  Step_size  Newton_decr.  Objective  log_det\n');

fGz = log(det(G'*diag(z)*G));
fz = -fGz - kappa*sum(log(z) + log(1-z));

fprintf('   0\t  -- \t     --   %10.3f  %10.3f\n', -fz, fGz);

for iter=1:MAXITER
    
    W = inv(G'*diag(z)*G);
    V = G*W*G';
%     V = G/(G'*diag(z)*G)*G';

    g = -2*real(diag(V)) - kappa*(1./z - 1./(1-z));
    H = 2*real(V.*V) + kappa*diag(1./(z.^2) + 1./((1-z).^2));

%     R = chol(H);
%     Hinvg = (R\(R'\g));
%     Hinv1 = (R\(R'\ones_M));
    Hinvg = H\g;
    Hinv1 = H\ones_M;
    dz = -Hinvg + ((ones_M'*Hinvg) / (ones_M'*Hinv1))*Hinv1;

    deczi = find(dz < 0);
    inczi = find(dz > 0);
    s = min([1; 0.99*[-z(deczi)./dz(deczi) ; (1-z(inczi))./dz(inczi)]]);

    while (1)
        zp = z + s*dz;
        fGzp = log(det(G'*diag(zp)*G));
        fzp = -fGzp - kappa*sum(log(zp) + log(1-zp));

        if (fzp <= fz + alpha*s*g'*dz)
            break;
        end
        s = beta*s;
    end
    z = zp; fz = fzp; fGz = fGzp;
    
    fprintf('%4d %10.3f %10.3f %10.3f %10.3f\n', iter, s, -g'*dz/2, -fz, fGz);

    if(-g'*dz/2 <= NT_TOL)
        break;
    end
end

zsort=sort(z); thres=zsort(M-K); zhat=(z>thres);
L = log(det(G'*diag(zhat)*G));
ztilde = z;
fGz = log(det(G'*diag(z)*G));
Utilde = fGz + 2*M*kappa;
fprintf('log_det: %f\n',L);
end
