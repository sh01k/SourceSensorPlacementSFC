function [idxc, idxr] = lmp_b_eim_n(U, L)
% Empirical interpolation method for joint source and sensor selection (broadband)
% INPUT
%     U                Functions to be approximated
%     L                Number of sources and sensors to be selected
% OUTPUT
%     idxc             Indexes of sensor locations
%     idxr             Indexes of source locations
% 
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

K = size(U,3);

idxr = zeros(L, 1);
idxc = zeros(L, 1);

Q = zeros(size(U, 1), L);
U = U ./ sqrt(repmat(sum(sum(abs(U).^2, 1),3), size(U, 1), 1, size(U, 3)));

for v = 1:L
    [k, ~] = max(abs(U),[],3);
    [m, idx1] = max(k);
    [~, idx2] = max(m);

    idxc(v) = idx2;
    idxr(v) = idx1(idx2);

    for kk=1:K
        Q(:, v, kk) = U(:, idx2, kk);

        U(:,:,kk) = U(:,:,kk) - Q(:, 1:v, kk) * (Q(idxr(1:v), 1:v, kk)\U(idxr(1:v), :, kk));
    end

end

end


