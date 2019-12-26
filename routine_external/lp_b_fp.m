function [select]= lp_b_fp(P,L)
% Source placement based on frame potential
% INPUT
%     P                Measurement matrix
%     L                Number of sources to be selected
% OUTPUT
%     select           Indexes of selected sources
%
% Modified from original code available here: 
%     https://github.com/jranieri/OptimalSensorPlacement
%
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

[N, M, J] = size(P);
select=1:M;

Pn = zeros(N, M, J);
for jj=1:J
    Pn(:,:,jj) = P(:,:,jj)*diag(1./sqrt(sum(abs(P(:,:,jj)).^2,1)));
end

%% Greedy algorithm

% First Step, find the best two columns
G = zeros(M,M);
for jj=1:J
    G = G + abs(Pn(:,:,jj)'*Pn(:,:,jj)).^2;
    G = G - eye(M);
end
Gsum=sum(G);

[~,ind] = max(Gsum);
Gsum=Gsum-G(select(ind),:);
select(ind)=[];

m=1;
% Find the column to eliminate one by one

while length(select)>L
    m=m+1;

    [~,ind] = max(Gsum(select));
    Gsum=Gsum-G(select(ind),:);
    select(ind)=[];
end

end
