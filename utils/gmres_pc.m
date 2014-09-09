clear

% load a random finite element matrix
nx = 4;
ny = 4;
A = gallery('wathen',nx,ny);

b = sum(A,2);


% set gmres parameters
tol = 1e-12;  
maxit = 1000;
restart = 30;


% solve with no preconditioner
%[x,flag,relres,iter] = gmres(A,b,restart,tol,maxit,[]);
%flag
%relres
%iter

% define a preconditioning matrix, matrix will invert this
B = diag(diag(A));

[x,flag,relres,iter] = gmres(A,b,restart,tol,maxit,B);
flag
relres
iter


% define a function to apply the inverse
Bmat = A; % Bmat NEED not be the same as A
pc_apply = @(r)preconditioner_invB(r,Bmat,b);

[x,flag,relres,iter] = gmres(A,b,restart,tol,maxit,pc_apply);
flag
relres
iter
