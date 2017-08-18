% N = 4; % dimension of the matrix

% A = complex(rand(N),rand(N)); % create random matrix

A = [1, 2, 3; 3, 1, 2; 2, 3, 1];
A = complex(A, zeros(3));

B = complex(eye(3), zeros(3));

[v1, e1] = eig(A, B);

[alpha, beta, v2] = mzggev(A, B);

e2 = alpha./beta;
