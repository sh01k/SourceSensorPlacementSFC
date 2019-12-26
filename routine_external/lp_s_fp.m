function [select]= lp_s_fp(P,L)
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
% Jun 2018 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

[N, M] = size(P);
select=1:M;

Pn = P*diag(1./sqrt(sum(abs(P).^2,1)));

%% Greedy algorithm

% First Step, find the best two columns
G = abs(Pn'*Pn).^2;
G = G - eye(M);
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
