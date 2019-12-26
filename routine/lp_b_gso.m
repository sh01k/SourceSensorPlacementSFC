function [select] = lp_b_gso(Z_cdt, K, d)
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

[M, L, J] = size(Z_cdt);
Z_cdt_idx = 1:L;

select = zeros(K,1);
Z = zeros(M,K,J);
V = zeros(M,K,J);
obj = zeros(K,J);
obj_av = zeros(K,1);

%Initialize
Fz = zeros(L,J);
for jj=1:J
    for ll=1:L
        pd = Z_cdt(:,ll,jj)'*d(:,jj)/(d(:,jj)'*d(:,jj))*d(:,jj);
        Fz(ll,jj) = norm(Z_cdt(:,ll,jj)-pd,2)/norm(Z_cdt(:,ll,jj),2);
    end
end

[~, select(1)] = max(sum(Fz,2)/J);
Z_cdt_idx(Z_cdt_idx == select(1))=[];
Z_cdt(:,Z_cdt_idx == select(1),:)=[];

Z(:,1,:) = Z_cdt(:,select(1),:);
for jj=1:J
    V(:,1,jj) = Z(:,1,jj)/norm(Z(:,1,jj),2);
    obj(1,jj) = norm(Z(:,1,jj),2);
end
obj_av(1) = sum(obj(1,:))/J;
fprintf('itr: 1, obj: %f, idx: %d\n', obj_av(1), select(1));

%Iterative selection
for kk=2:K
    r = zeros(M, K-kk, J);
    r_norm = zeros(K-kk, J);
    for jj=1:J
        for ll=1:L-kk
            p = (repmat(V(:,1:(kk-1),jj)'*Z_cdt(:,ll,jj),1,M).').*V(:,1:(kk-1),jj);
            r(:,ll,jj) = Z_cdt(:,ll,jj) - sum(p,2);
            r_norm(ll,jj) = norm(r(:,ll,jj),2);
        end
    end
    
    [~,idx] = max(sum(r_norm,2)/J);
    select(kk) = Z_cdt_idx(idx);
    
    Z(:,kk,:) = Z_cdt(:,select(kk),:);
    for jj=1:J
        V(:,kk,jj) = r(:,idx,jj)/norm(r(:,idx,jj),2);
        obj(kk,jj) = r_norm(idx,jj);
    end
    obj_av(kk) = sum(obj(kk,:),2)/J;
    
    Z_cdt_idx(Z_cdt_idx==select(kk))=[];
    Z_cdt(:,Z_cdt_idx==select(kk),:)=[];
    
    fprintf('itr: %d, obj: %f, idx: %d\n',kk, obj_av(kk), select(kk));
end

end