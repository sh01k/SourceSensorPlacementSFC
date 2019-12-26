function [z_loc, L_loc] = mp_s_det_loc(G, zhat)
% Local optimization method for sensor selection
%     http://ieeexplore.ieee.org/document/4663892/
% INPUT
%     G                Measurement matrix
%     zhat             Initial z of binary
% OUTPUT
%     z_loc            Estimate z
%     L_loc            Objective value
% 
% Modified from original code available here: 
%     www.stanford.edu/~boyd/papers/sensor_selection.html
% 
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

[M, N] = size(G);
K = sum(zhat);
zn = zhat;
iteration = 0; 
N_loc = 10*M^3/N^2; 
swapstaken = 0;

hatSig = inv(G'*diag(zn)*G);

while(iteration < N_loc)
    flag = 0;
    S = find(zn == 1); nS = find(zn == 0);   
    
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
                zn(out) = 0; zn(in) = 1;
                flag = 1;
                break;
            end
        end
        if( flag == 1), break; end
    end
    if(flag==0), break; end
end

fprintf('\nLocal optimization:\n');
if(flag==1)
   fprintf('Maximum Iteration = %d reached. Terminated.\n', N_loc);
end
fprintf('Swaps checked = %d, swaps taken = %d, [K(M-K) = %d]\n', iteration, swapstaken, length(S)*length(nS));
z_loc = zn;
L_loc = log(det(G'*diag(z_loc)*G));

end
