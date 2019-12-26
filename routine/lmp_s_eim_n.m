function [idxc, idxr] = lmp_s_eim_n(U, L)
% Empirical interpolation method for joint source and sensor selection (single frequency)
% INPUT
%     U                Functions to be approximated
%     L                Number of sources and sensors to be selected
% OUTPUT
%     idxc             Indexes of sensor locations
%     idxr             Indexes of source locations
% 
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

idxr = zeros(L, 1);
idxc = zeros(L, 1);

Q = zeros(size(U, 1), L);
U = U ./ sqrt(repmat(sum(abs(U).^2, 1), size(U, 1), 1));

for v = 1:L
    [m, idx1] = max(abs(U));
    [~, idx2] = max(m);

    idxc(v) = idx2;
    idxr(v) = idx1(idx2);
    Q(:, v) = U(:, idx2);

    U = U - Q(:, 1:v) * (Q(idxr(1:v), 1:v)\U(idxr(1:v), :));
end

end


