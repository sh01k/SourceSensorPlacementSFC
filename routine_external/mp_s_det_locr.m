function [z_loc, L_loc] = mp_s_det_locr(G, K, zast, threshold)
% Local optimization method for sensor selection
%     http://ieeexplore.ieee.org/document/4663892/
% INPUT
%     G                Measurement matrix
%     K                Number of sensors to be selected
%     zast             Initial z of real value
%     thresold         Constant parameter for thresholding
% OUTPUT
%     z_loc            Estimate z
%     L_loc            Objective value
% 
% Modified from original code available here: 
%     www.stanford.edu/~boyd/papers/sensor_selection.html
% 
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

if(nargin < 4); threshold = 0.5; end

[M, N] = size(G);
N_loc = 10*M^3/N^2;

[z_s, idx_s] = sort(zast); 
thres = z_s(M-K); zhat=(zast>thres); zhat_s = (z_s > thres);
toswap01 = (abs(z_s - 0.5) <= threshold);
S = idx_s(toswap01 & zhat_s); 
nS = flipud(idx_s(toswap01 & ~zhat_s)); 

fprintf('\nLocal optimization (ordering and thresholding):\n');
fprintf('Number of sensors to chosen to swap: %d ', sum(toswap01));
fprintf('selected: %d, not selected %d\n', length(S), length(nS));

iteration = 0; 
swapstaken = 0;

hatSig = inv(G'*diag(zhat)*G);

while(iteration < N_loc)
    flag = 0;
    
    for out = S' %matlab does not allow to say for out = S
        for in = nS' 
            iteration = iteration + 1;
            v1 = [G(out,:); G(in, :)];
            v2 = [-G(out,:); G(in, :)]';
            S_2by2 = eye(2) + v1*hatSig*v2;
            volchangefactor = det(S_2by2);
            if(volchangefactor > 1)
                fprintf('Value= %f\t Swap: OUT %d\tIN %d\n', volchangefactor, out, in);
                swapstaken = swapstaken + 1;
                hatSig = hatSig - hatSig*v2*inv(S_2by2)*v1*hatSig;
                zhat_s(idx_s == out) = 0; zhat_s(idx_s == in) = 1;
                S = idx_s(toswap01 & zhat_s); 
                nS = flipud(idx_s(toswap01 & ~zhat_s)); 
                flag = 1;
                break;
            end
        end
        if( flag == 1), break; end
    end
    if(flag==0), break; end
end


if(flag==1)
   fprintf('Maximum Iteration = %d reached. Terminated.\n', N_loc);
end
fprintf('Swaps checked = %d, swaps taken = %d, len(S)*len(nS) = %d]\n', iteration, swapstaken, length(S)*length(nS));
z_loc = zeros(M,1); z_loc(idx_s(zhat_s)) = 1;
L_loc = log(det(G'*diag(z_loc)*G));
fprintf('log_det: %f\n',L_loc);

z_loc = logical(z_loc);
end
