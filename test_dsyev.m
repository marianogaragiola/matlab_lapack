N = 4; % dimension of the matrix

A = rand(N); % create random matrix

A = A + A'; % simetrization of the matrix

e1 = eig(A); % calculation of the eigenvalue 

e2 = mdsyev(A);
