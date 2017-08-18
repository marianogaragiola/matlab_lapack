N = 4; % dimension of the matrix

A = rand(N); % create random matrix

A = A + A'; % simetrization of the matrix

B = eye(4);

[v1, e1] = eig(A); % calculation of the eigenvalue

[e2, v2] = mdsygvx(A, B, 0, 0);
