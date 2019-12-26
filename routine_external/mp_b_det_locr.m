function [z_loc, L_loc] = mp_b_det_locr(G, K, zast, threshold)
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

[M, N, J] = size(G);
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

hatSig = zeros(K,K,J);
dethatSig = zeros(J,1);
for jj=1:J
%     hatSig(:,:,jj) = inv(G(:,:,jj)'*diag(zhat)*G(:,:,jj));
    [U_Sig,S_Sig,V_Sig]=svd(G(:,:,jj)'*diag(zhat)*G(:,:,jj));
    S_Sig(S_Sig~=0)=S_Sig(S_Sig~=0).^-1;
    hatSig(:,:,jj)=V_Sig*S_Sig*U_Sig';
    
    dethatSig(jj) = det(hatSig(:,:,jj));
end

while(iteration < N_loc)
    flag = 0;
    
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
L_loc = 0;
for jj=1:J
    L_loc = L_loc + log(det(G(:,:,jj)'*diag(z_loc)*G(:,:,jj)));
end
L_loc = L_loc/J;

z_loc = logical(z_loc);
fprintf('log_det: %f\n',L_loc);
end


