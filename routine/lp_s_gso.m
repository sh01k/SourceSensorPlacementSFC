function [select] = lp_s_gso(Z_cdt, K, d)
% Source placement by Gram-Schmidt orthogonalization
%     http://ieeexplore.ieee.org/document/748126/
% INPUT
%     Z_cdt            Transfer function matrix of candidate locations
%     K                Number of sources to be selected
%     d                Desired sound field for initialization
% OUTPUT
%     select           Indexes of source locations
% 
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

[M, L] = size(Z_cdt);
Z_cdt_idx = 1:L;

select = zeros(K,1);
Z = zeros(M,K);
V = zeros(M,K);
obj = zeros(K,1);

%Initialize
Fz = zeros(L,1);
for ll=1:L
    pd = Z_cdt(:,ll)'*d/(d'*d)*d;
    Fz(ll) = norm(Z_cdt(:,ll)-pd,2)/norm(Z_cdt(:,ll),2);
end

[~, select(1)] = max(Fz);
Z_cdt_idx(Z_cdt_idx == select(1))=[];
Z_cdt(:,Z_cdt_idx == select(1))=[];

Z(:,1) = Z_cdt(:,select(1));
V(:,1) = Z(:,1)/norm(Z(:,1),2);
obj(1) = norm(Z(:,1),2);
fprintf('itr: 1, obj: %f, idx: %d\n', obj(1), select(1));

%Iterative selection
for kk=2:K
    r = zeros(M, K-kk);
    r_norm = zeros(K-kk, 1);
    for ll=1:L-kk
        p = (repmat(V(:,1:(kk-1))'*Z_cdt(:,ll),1,M).').*V(:,1:(kk-1));
        r(:,ll) = Z_cdt(:,ll) - sum(p,2);
        r_norm(ll) = norm(r(:,ll),2);
    end
    
    [~,idx] = max(r_norm);
    select(kk) = Z_cdt_idx(idx);
    
    Z(:,kk) = Z_cdt(:,select(kk));
    V(:,kk) = r(:,idx)/norm(r(:,idx),2);
    obj(kk) = r_norm(idx);
    
    Z_cdt_idx(Z_cdt_idx==select(kk))=[];
    Z_cdt(:,Z_cdt_idx==select(kk))=[];
    
    fprintf('itr: %d, obj: %f, idx: %d\n',kk, obj(kk), select(kk));
end

end