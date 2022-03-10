function err = Poi2D(k)
tic
N = 2^k;
h = 1/N;

% Functions of vector field and analytical solutions
f=@(x,y) -(sin(pi*x)).*(sin(pi*y));

ana=@(x,y) (sin(pi*x).*sin(pi*y))/(-2*pi*pi);

% Prepare Grid
[X,Y] = meshgrid(h:h:1-h ,h:h:1-h);

% Construct Matrices
FD=(1/h^2)*spdiags([-ones(N-1,1) 2*ones(N-1,1) -ones(N-1,1)],-1:1,N-1,N-1);
ID = sparse(eye(N-1));

% Kron-Products
A = kron(FD,ID) + kron(ID,FD);
b = reshape(f(X,Y),[(N-1)^2,1]);

% Solving the System  
dsol = (A\b) ;

% Matrix of Solution on Grid
grid_sol = reshape(dsol, [N-1, N-1]);

% Compute Errors
ana_sol = reshape(ana(X,Y), [(N-1)^2, 1]);

err = norm(dsol-ana_sol)/sqrt(length(dsol));

surf(X,Y, grid_sol)
toc
end