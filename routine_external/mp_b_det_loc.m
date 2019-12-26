function [z_loc, L_loc] = mp_b_det_loc(G, zhat)
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

[M, N, J] = size(G);
K = sum(zhat);
zn = zhat;
iteration = 0; 
N_loc = 10*M^3/N^2; 
swapstaken = 0;

hatSig = zeros(K,K,J);
dethatSig = zeros(J,1);
for jj=1:J
    hatSig(:,:,jj) = inv(G(:,:,jj)'*diag(zn)*G(:,:,jj));
    dethatSig(jj) = det(hatSig(:,:,jj));
end

while(iteration < N_loc)
    flag = 0;
    S = find(zn == 1); nS = find(zn == 0);   
    
    for out = S' %matlab does not allow to say for out = S
        for in = nS' 
            iteration = iteration + 1;
            v1 = zeros(2, K, J);
            v2 = zeros(K, 2, J);
            S_2by2 = zeros(2,2,J);
            detS_2by2 = zeros(J,1);
            for jj=1:J
                v1(:,:,jj) = [G(out,:,jj); G(in, :,jj)];
                v2(:,:,jj) = [-G(out,:,jj); G(in, :,jj)]';
                S_2by2(:,:,jj) = eye(2) + v1(:,:,jj)*hatSig(:,:,jj)*v2(:,:,jj);
                detS_2by2(jj) = det(S_2by2(:,:,jj));
            end
            volchangefactor = sum(dethatSig.*detS_2by2)/sum(dethatSig);
            if(volchangefactor > 1)
                fprintf('Value= %f\t Swap: OUT %d\tIN %d\n', volchangefactor, out, in);
                swapstaken = swapstaken + 1;
                for jj=1:J
                    hatSig(:,:,jj) = hatSig(:,:,jj) - hatSig(:,:,jj)*v2(:,:,jj)*inv(S_2by2(:,:,jj))*v1(:,:,jj)*hatSig(:,:,jj);
                    dethatSig(jj) = det(hatSig(:,:,jj));
                end
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
L_loc = 0;
for jj=1:J
    L_loc = L_loc + log(det(G(:,:,jj)'*diag(z_loc)*G(:,:,jj)));
end
L_loc = L_loc/J;

end
