function binconn = binarize_conn(conn)
% 
% This function binarize weighted connectivity matrix with threshold of 95%
% and return binarized connectivity matrix.
% 
% Inputs :  conn    - Weighted connectivity matrix  (n x n)
% Outputs:  binconn - Binarized connectivity matrix (n x n)
% 

th = 5;
Nroi = size(conn, 1);
binconn = zeros(Nroi, Nroi);

W = conn;
if isequal(W,W.')	% if symmetric matrix
    W=triu(W);      % ensure symmetry is preserved
    ud=2;           % halve number of removed links
else
    ud=1;
end
    
ind=find(W);                            % find all links
E=sortrows([ind W(ind)], -2);           % sort by magnitude
en=round((Nroi^2-Nroi)*0.01*th/ud);     % number of links to be preserved. the value in front of '/ud' will be changed
W(E(en+1:end,1))=0;                     % apply threshold
    
% if symmetric, reconstruct symmetry
if ud==2
    W=W+W.';
end

% binarize
W=double(W~=0);
binconn = W;
end
