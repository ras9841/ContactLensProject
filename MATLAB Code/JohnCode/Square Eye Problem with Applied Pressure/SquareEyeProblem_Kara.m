%% Saran Wrap Problem Statement
% Solving the differential equation/BVP:
% a<=x<=b;
% c<=z<=d;
% 
% Left side: u(a,z) = 0; w(a,z) = 0; (Dirichlet conditions)
% Right side: u(b,z) = constant; w(b,z) = 0; (Dirichlet conditions)
% Top side: (sigma)*u_x(x,d) + w_z(x,d) = 0 (Result of normal stress being
% zero in y-direction)
% Top side: u_z(x,d) + w_x(x,d) = 0 (Result of shear stress in sigma_{zx} = 0
% Bottom side: w(x,c) = 0 (Dirichlet condition)
% Bottom side: u_z(x,c) = 0 (Neumann condition)
% PDE-1: (u_xx + u_zz) + (1/(1-sigma))*(u_xx + w_xz) = 0;
% PDE-2: (w_xx + w_zz) + (1/(1-sigma))*(u_xz + w_zz) = 0;

% Notice: number of grid points in the x-direction and the z-direction are
% the same (n) for both u(x,z) and w(x,z).
%% Creating the derivative matrix to approximate u_zz
% remember u(x,z) = 0 for all boundary points
% this just means the vec_b is all zeros

tic;

% Poisson's Ratio (must be on interval of [-1,1/2))
sigma = 0.49; % sigma = 0.4998 works well

n = 5; % number of POINTS in x-direction and z-direction for both 
% variables u and w.

% Young's Modulus
E = 1;

k1 = 1/(1-sigma);
k2 = 1/(1-2*sigma);
k = k2;

% x is on interval a<x<b
% z is on interval c<z<d
a = 0; b = 1;
c = 0; d = 1;

hx = (b-a)/(n-1); % step size in the x-direction: 
% there are n points to be evaluated in the x-direction
hz = (d-c)/(n-1); % step size in the z-direction:
% there are n points to be evaluated in the z-direction

X = linspace(a,b,n);
Z = linspace(c,d,n);
[X,Z]=meshgrid(X,Z);

e = ones(n,1); % a column vector with n ones
D2z_u = sparse(1:n, [2:n 1], e, n, n); % sparse matrix with ones down upper-
% off-diagonal and a one in the lower left corner
D2z_u = D2z_u+D2z_u'; % D2 now has ones down both upper- and lower-off-diagonals
% and ones in the lower left and upper right corners
D2z_u(1:n+1:n^2) = -2; % steps through each row and changes diagonal elements to 2
D2z_u = D2z_u/(hz^2); % divided each entry by the step size squared: (h_z)^2
% at this point the entire matrix is correct except with respect to the 
% boundary conditions.
D2z_u(1,:) = [1 zeros(1, n-1)]; % changes the first row to a one followed by all zeros
D2z_u(end, :) = [zeros(1, n-1) 1]; % changes the last row to all zeros followed by a one

I_n = eye(n); % nxn identity matrix
u_zz = kron(I_n, D2z_u); % using kronecker product to create matrix that almost approximates u_zz
% u_zz is an (n^2)x(n^2) matrix
u_zz(1:n, 1:n) = eye(n); % changes the upper-left block from D2z to an identity matrix
u_zz(end-n+1:end, end-n+1:end) = eye(n); % changes the lower-right block from D2z to an identity matrix
% now u_zz finally approximates the second derivative u_zz

%% Approximating D2x_u with a matrix
% approach is generally the same as for u_zz above
e = ones(n,1); % a column vector with n ones
D2x_u = sparse(1:n, [2:n 1], e, n, n); % sparse matrix with ones down upper-
% off-diagonal and a one in the lower left corner
D2x_u = D2x_u+D2x_u'; % D2 now has ones down both upper- and lower-off-diagonals
% and ones in the lower left and upper right corners
D2x_u(1:n+1:n^2) = -2; % steps through each row and changes diagonal elements to 2
D2x_u = ((k+1)/hx^2)*D2x_u; % Multiply by the constant: (k+1)/(h_x)^2
% at this point the entire matrix is correct except with respect to the 
% boundary conditions.
D2x_u(1,:) = [1 zeros(1, n-1)]; % changes the first row to a one followed by all zeros
D2x_u(end, :) = [zeros(1, n-1) 1]; % changes the last row to all zeros followed by a one

% I_n = eye(n); % nxn identity matrix
u_xx = kron(D2x_u, I_n); % using kronecker product to create matrix that almost approximates u_xx
% u_zz is an (n^2)x(n^2) matrix
u_xx(1:n, 1:n) = eye(n); % changes the upper-left block from D2x to an identity matrix
u_xx(end-n+1:end, end-n+1:end) = eye(n); % changes the lower-right block from D2x to an identity matrix
% this is not the same as u_zz because u_xx is a tri-diagonal block matrix,
% not just a tri-diagonal matrix like u_zz
% we also need to change certain rows to respect boundary conditions
% such rows are all zeros except a one on the diagonal
% It is straightforward to replace each row's 3 appropriate elements using
% a for loop. 0
for i=1:(n-2),
    upperRow = i*n + 1;
    u_xx(upperRow, (i-1)*n+1)=0;
    u_xx(upperRow, (i)*n+1)=0; % puts a zero on the diagonal of u_xx
    u_xx(upperRow, (i+1)*n+1)=0;
    
    lowerRow = i*n + n;
    u_xx(lowerRow, (i-1)*n+n)=0;
    u_xx(lowerRow, (i)*n+n)=0; % puts a zero on the diagonal of u_xx
    u_xx(lowerRow, (i+1)*n+n)=0;
end;
% now u_xx finally approximates the second derivative u_xx
%%
%%
%% Creating the derivative matrix to approximate w_zz
% remember u(x,z) = 0 for all boundary points EXCEPT top where w(x,d) =
% f(x) is a displacement function

% x is on interval a<x<b
% z is on interval c<z<d


D2z_w = sparse(1:n, [2:n 1], e, n, n); % sparse matrix with ones down upper-
% off-diagonal and a one in the lower left corner
D2z_w = D2z_w+D2z_w'; % D2 now has ones down both upper- and lower-off-diagonals
% and ones in the lower left and upper right corners
D2z_w(1:n+1:n^2) = -2; % steps through each row and changes diagonal elements to 2
D2z_w = ((k+1)/hz^2)*D2z_w; % multiply by the constant (k+1)/(h_z)^2
% at this point the entire matrix is correct except with respect to the 
% boundary conditions.
D2z_w(1,:) = [1 zeros(1, n-1)]; % changes the first row to a one followed by all zeros
D2z_w(end, :) = [zeros(1, n-1) 1]; % changes the last row to all zeros followed by a one

I_n = eye(n); % nxn identity matrix
w_zz = kron(I_n, D2z_w); % using kronecker product to create matrix that almost approximates u_zz
% u_zz is an (n^2)x(n^2) matrix
w_zz(1:n, 1:n) = eye(n); % changes the upper-left block from D2z to an identity matrix
w_zz(end-n+1:end, end-n+1:end) = eye(n); % changes the lower-right block from D2z to an identity matrix
% now u_zz finally approximates the second derivative u_zz

%% Approximating D2x_w with a matrix
% approach is generally the same as for u_zz above
e = ones(n,1); % a column vector with n ones
D2x_w = sparse(1:n, [2:n 1], e, n, n); % sparse matrix with ones down upper-
% off-diagonal and a one in the lower left corner
D2x_w = D2x_w+D2x_w'; % D2 now has ones down both upper- and lower-off-diagonals
% and ones in the lower left and upper right corners
D2x_w(1:n+1:n^2) = -2; % steps through each row and changes diagonal elements to 2
% there are (n_x + 1) points to be evaluated in the x-direction
D2x_w = D2x_w/(hx^2); % divided each entry by the step size squared: (h_x)^2
% at this point the entire matrix is correct except with respect to the 
% boundary conditions.
D2x_w(1,:) = [1 zeros(1, n-1)]; % changes the first row to a one followed by all zeros
D2x_w(end, :) = [zeros(1, n-1) 1]; % changes the last row to all zeros followed by a one

% I_n = eye(n); % nxn identity matrix
w_xx = kron(D2x_w, I_n); % using kronecker product to create matrix that almost approximates u_xx
% u_zz is an (n^2)x(n^2) matrix
w_xx(1:n, 1:n) = eye(n); % changes the upper-left block from D2x to an identity matrix
w_xx(end-n+1:end, end-n+1:end) = eye(n); % changes the lower-right block from D2x to an identity matrix
% this is not the same as u_zz because u_xx is a tri-diagonal block matrix,
% not just a tri-diagonal matrix like u_zz
% we also need to change certain rows to respect boundary conditions
% such rows are all zeros except a one on the diagonal
% It is straightforward to replace each row's 3 appropriate elements using
% a for loop.
for i=1:(n-2),
    upperRow = i*n + 1;
    w_xx(upperRow, (i-1)*n+1)=0;
    w_xx(upperRow, (i)*n+1)=0; % puts a zero on the diagonal of w_xx
    w_xx(upperRow, (i+1)*n+1)=0;
    
    lowerRow = i*n + n;
    w_xx(lowerRow, (i-1)*n+n)=0;
    w_xx(lowerRow, (i)*n+n)=0; % puts a zero on the diagonal of w_xx
    w_xx(lowerRow, (i+1)*n+n)=0;
end;
% now u_xx finally approximates the second derivative u_xx
%% Setting up the Laplacians
Lap_u = u_xx + u_zz; 
% boundary conditions have been taken care of earlier but L has a sum of
% Identity matrices in the blocks of the upper-left and lower-right
% corners. They can simply be halved or replaced with new matrices.
% also, the constant in front of u_xx was taken care of earlier.
Lap_u(1:n, 1:n) = eye(n);
Lap_u(end-n+1:end, end-n+1:end) = eye(n);
% Also, the first and last entries of the nxn blocks on the diagonals need
% to be halved. 



Lap_w = w_xx + w_zz; 
% boundary conditions have been taken care of earlier but L has a sum of
% Identity matrices in the blocks of the upper-left and lower-right
% corners. They can simply be halved or replaced with new matrices.
% also, the constant in front of u_xx was taken care of earlier.
Lap_w(1:n, 1:n) = eye(n);
Lap_w(end-n+1:end, end-n+1:end) = eye(n);

%% Coupling the PDEs - Approximating D2xz with a matrix
% Instead of making a block diagonal matrix with only Lap_u and Lap_w, we
% need to couple the system by having terms in the upper right and lower
% left quadrants of what would be a large (2n^2)x(2n^2) matrix. The
% creation is very similar to generating D2x_u and D2x_w. 

e = ones(n,1); % a column vector with n ones
m1 = sparse( 1:n, [2:n 1], e, n, n ); % sparse matrix with ones down the upper off-diagonal and a one in the lower left corner
m1 = m1'; % m1 now has only ones down its lower off-diagonal and a one in the upper right corner.
m2 = m1-m1'; % m2 has +1 down its lower off-diagonal and -1 down its upper off diagonal. Need to clear the bottom rows.
m2( 1, : ) = zeros( 1, n ); % clears m2's top row
m2( n, : ) = zeros( 1, n ); % clears m2's bottom row
m3 = -m2; % m3 has -1 down its lower off-diagonal and +1 down its upper off diagonal, except the top and bottom rows
m4 = horzcat( horzcat( m2, zeros( n, n )), m3 );
m5 = spalloc(n^2,n^2,2*(n-2)); % each of the large "coupling matrices" will be (n^2)x(n^2) with exactly 2*(n-2) nonzero entries
for i=2:(n-1),
    m5( (i-1)*n+1:i*n, (i-2)*n+1:(i+1)*n ) = m4; %#ok<SPRIX>
end;
% m5 now has m2 matrices down its lower off-diagonal blocks and m3 matrices
% down its upper off-diagonal blocks. 
% The top and bottom row blocks, which you can think of as the nx(n^2) 
% block on top and bottom of the matrix m5, are all zeros.
% Need to multiply all entries by the constant from the PDE
D2xz = (k/(4*hx*hz))*m5; % D2xz now approximates the coupling part of the PDEs. It works for both the first and second PDE

%% Solving the Laplacian (and displaying it?)
% to solve, set up as (Lap_u + kdD/dx : Lap_w +kdD/dz)(u .. w) = vec_b 
%Lap_uw = blkdiag( Lap_u, Lap_w );
FirstPDE = horzcat( Lap_u, D2xz );
SecondPDE = horzcat( D2xz, Lap_w );
A = vertcat( FirstPDE, SecondPDE );


%% %% Correcting the matrix A for appropriate BCs
%% Fixing A for top conditions 
% \sigma_{zz} = p(x)  ==>  (E*sigma)/((1+sigma)(1-2*sigma))*( (1-sigma)w_z + sigma(u_x)) = p(x)
% ==> t1*w_z + t2*u_x = p(x) where t1 and t2 are corresponding contstants above;
% 
% Also, \sigma_{xz} = 0  ==>  u_x + w_z = 0;
% The boundary conditions do not include the values at the corners, that is
% something to take care of for the roots of the function p(x)


r1 = spalloc( 1, 3*n+n^2, 7 );
r2 = spalloc( 1, 3*n+n^2, 7 );
r3 = spalloc( 1, 3*n+n^2, 7 );
% Note: r1, r2, r3 are instantiated as being row vectors that are
% 1-by-(3*n+n^2) with memory allocated for at most 7 nonzero entries because
% they are sparse matrices. However, none of them will actually have 7
% entries (r1 has to do with u_xx and will have 3 nonzero entries, r2 has
% to do with u_z and  will have 2 non-zero entries, and r3 has to do with
% w_x and will have 2 nonzero entries). In total though, they are to be
% added and would have at most 7 entries in their sum, but here the u_i,j
% (u_center) term from u_z and u_xx will add together for 6 nonzero terms.
% The 7th memory spot is just there to ensure no memory conflicts will
% arise.


% it can be shown that when V=V_y=V_xy=V_yz=0,
% w_xz = t1*u_xx + t2*p'(x) where p'(x) is the derivative of the pressure
% function p(x). To show this, differentiate the sigma_zz normal pressure
% with respect to x and solve for w_xz.

if k == k1
    alpha1 = ( -2*E ) / ( (1+sigma)*hx^2 );
    alpha2 = ( -2*E ) / ( (1+sigma)*hz^2 );
    alpha3 = ( E ) / ( (1+sigma)*hx*hz );
elseif k == k2
    alpha1 = ( -E*(2-sigma) ) / ( (1+sigma)*hx^2 );
    alpha2 = ( -2*E*(1-sigma) ) / ( (1+sigma)*hz^2 );
    alpha3 = ( E*(1-sigma) ) / ( (1+sigma)*hx*hz );
end;

r1( 1, n ) = 1;
r1( 1, 2*n ) = -2;
r1( 1, 3*n ) = 1;
R1 = alpha1*r1;

r2( 1, (2*n-1):(2*n)) = [1 -1];
R2 = alpha2*r2;

r3( 1, n+n^2 ) = -1;
r3( 1, 3*n+n^2 ) = 1;
R3 = alpha3*r3;

% "tbcrvu" stands for "top boundary conditions row vector [from the
% perspective of] u"
tbcrvu = R1+R2+R3;
for i = 2:(n-1),
    A( i*n, :) = zeros(1,2*n^2);
    A( i*n, ((i-2)*n+1):((i-2)*n+3*n+n^2)) = tbcrvu;
    % this takes care of the top boundary conditions from U's perspective
end;

%----------------------------------------------------------------------%

r4 = spalloc( 1, 3*n+n^2, 7 );
r5 = spalloc( 1, 3*n+n^2, 7 );
r6 = spalloc( 1, 3*n+n^2, 7 );
% See the comments for r1, r2, and r3 to understand why the maximum number
% of nonzero entries is set to 7 for each of r4, r5, r6


% it can be shown that when V=V_y=V_xy=V_yz=0,
% w_z = t3*u_x + t4*p(x) where p(x) is the pressure function. To do this,
% express sigma_zz (the normal stress in the z-direction) in terms of w_z, 
% u_x, and p(x), then solve for w_z.

if k == k1
    alpha4 = ( E*sigma*hz ) / (2*(1-sigma^2)*(2-sigma)*hx^2 );
    alpha5 = ( -E ) / ( (1-sigma^2)*hz );
    alpha6 = ( E*sigma ) / ( 2*(1-sigma^2)*hx );
elseif k == k2
    alpha4 = ( E*sigma*hz ) / ( 2*(1+sigma)*(1-2*sigma)*hx^2 );
    alpha5 = ( -E*(1-sigma) ) / ( (1+sigma)*(1-2*sigma)*hz );
    alpha6 = ( E*sigma ) / ( 2*(1+sigma)*(1-2*sigma)*hx );
end;

r4( 1, n+n^2 ) = 1;
r4( 1, 2*n+n^2 ) = -2;
r4( 1, 3*n+n^2 ) = 1;
R4 = alpha4*r4;

r5( 1, (2*n-1+n^2):(2*n+n^2) ) = [1 -1];
R5 = alpha5*r5;

r6( 1, n ) = -1;
r6( 1, 3*n ) = 1;
R6 = alpha6*r6;

% tbcrvw stands for "top boundary conditions row vector [from the
% perspective of] w"
tbcrvw = R4+R5+R6;
for i = 2:(n-1),
    A( i*n+n^2, : ) = zeros(1,2*n^2);
    A( i*n+n^2, ((i-2)*n+1):(n^2+3*n+(i-2)*n) ) = tbcrvw;
    % this takes care of the top boundary conditions from W's perspective
end;
%% Boundary Conditions (BCs) vector
% Boundary conditions for the corners (where displacement functions on the 
% sides supercede the derivative conditions on the top)


vec_p = zeros( n, 1 ); % vec_p == \vec{p}(x) which is the pressure function vector
for i = 1:n,
    vec_p( i, 1 ) = SquareEyeStressFunction( a+(i-1)*hx );
end;
vec_p_prime = zeros( n, 1 ); % p_vec_prime == \vec{p}'(x) which is the 
% derivative of the presssure function vector with respect to x
vec_p_prime( 1, 1 ) = 2*( vec_p(2,1) - vec_p(1,1) ); % top left corner is forward/backward difference
vec_p_prime( n, 1 ) = 2*( vec_p(n,1) - vec_p(n-1,1) ); % top right corner is forward/backward difference
for i = 2:(n-1),
    vec_p_prime( i, 1 ) = vec_p( i+1, 1 ) - vec_p( i-1, 1 );
end;
vec_p_prime = (1/(2*hx))*vec_p_prime;

%vec_b_u = zeros(n,n);
%vec_b_u = -(k+1)*((2*pi)^2)*cos(2*pi*X);
%vec_b_u = vec_b_u(:);
for i=1:n,
    vec_b_u( n*i, 1 ) = vec_p_prime( i, 1 ); % "top" (z=d); exact value of the pressure function
   % vec_b_u((i-1)*n+1, 1) = cos(2*pi*(a+(i-1)*hx));
  %  vec_b_u(i, 1) = 1;
 %   vec_b_u((n-1)*n+i, 1) = 1;
%     vec_b_u((i-1)*n+1, 1) = SquareEyeDisplacementFunction_U('bottom', a + (i-1)*hx, c);
%     vec_b_u(i, 1) = SquareEyeDisplacementFunction_U('left', a, c + (i-1)*hz);
%     vec_b_u((n-1)*n+i, 1) = SquareEyeDisplacementFunction_U('right', b, c + (i-1)*hz);
end;
%vec_b_u( n, 1 ) = 0; % causes top left corner to be zero (since U=W=0 on the left)
%vec_b_u( n^2, 1 ) = 0; % cause top right corner to be zero (since U=W=0 on the right)
    
    
%vec_b_w = zeros(n,n);
%vec_b_w = -(k+1)*((2*pi)^2)*cos(2*pi*Z);
%vec_b_w = vec_b_w(:);

for i=1:n,
    vec_b_w(n*i, 1) = vec_p( i, 1 ); % "top" (z=d); approximation of the derivative of the pressure function
  %  vec_b_w((i-1)*n+1, 1) = 1;
 %   vec_b_w(i, 1) = cos(2*pi*(c + (i-1)*hz));
%    vec_b_w((n-1)*n+i, 1) = cos(2*pi*(c + (i-1)*hz));
    %     vec_b_w((i-1)*n+1, 1) = SquareEyeDisplacementFunction_W('bottom', a + (i-1)*hx, c);
%     vec_b_w(i, 1) = SquareEyeDisplacementFunction_W( 'left', a, c + (i-1)*hz);
%     vec_b_w((n-1)*n+i, 1) = SaranWrapDisplacementFunction_W_only('right', b, c + (i-1)*hz);
end;
vec_b_w( n, 1 ) = 0; % causes top left corner to be zero (since U=W=0 on the left)
vec_b_w( n^2, 1 ) = 0; % cause top right corner to be zero (since U=W=0 on the right)
vec_b_uw = vertcat( vec_b_u, vec_b_w );


%% Solving the problem
uw1 = A\vec_b_uw; % solution is found here
uw2 = reshape( uw1, n, 2*n );
u_sol = uw2( 1:n, 1:n );
w_sol = uw2( 1:n, (n+1):2*n );
%w_sol = flipud( w_sol );

% %%% Iterative approach
% A_diag = diag(diag(A));
% A_diag_slash_b = A_diag\vec_b_uw;
% A_diag_slash_A = A_diag\(A);
% 
% x_iter = zeros(2*n*n,1);
% 
% for i=1:100
%    x_iter = x_iter + 0*(A_diag_slash_b - A_diag_slash_A*x_iter);   
% end
% 
% uw2 = reshape( x_iter, n, 2*n );
% u_sol = uw2( 1:n, 1:n );
% w_sol = uw2( 1:n, (n+1):2*n );



figure(1);
surf( X,Z, u_sol ) % plot of u's solution
%shading interp;
title( 'U(x,z)', 'FontSize', 48,...
       'FontWeight','bold');
xlabel('X axis','FontSize', 36,...
       'FontWeight','bold');
ylabel( 'Z axis','FontSize', 36,...
       'FontWeight','bold' );
zlabel( 'Displacement in the X-direction', 'FontSize', 18,...
       'FontWeight','bold' );


   
   
   
figure(2);
surf( X, Z, w_sol ) % plot of w's solution
title( 'W(x,z)', 'FontSize', 48,...
       'FontWeight','bold');
xlabel('X axis','FontSize', 36,...
       'FontWeight','bold');
ylabel( 'Z axis','FontSize', 36,...
       'FontWeight','bold' );
zlabel( 'Displacement in the Z-direction', 'FontSize', 18,...
       'FontWeight','bold' );
%shading interp;




figure(3);
surf( X+u_sol, zeros(n,n), Z+w_sol )
%shading interp;
title( 'Combined Solution', 'FontSize', 48,...
       'FontWeight','bold');
xlabel('X axis','FontSize', 36,...
       'FontWeight','bold');
ylabel( 'Y axis','FontSize', 36,...
       'FontWeight','bold' );
zlabel( 'Z axis', 'FontSize', 36,...
       'FontWeight','bold' );

% % % s1 = ((1-2*sigma^2)/(1-sigma))*(1/(2*hx));
% % % s2 = ((sigma*(1-2*sigma))/(1-sigma))*(1/(2*hz));
% % % s3 = 1/(2*hz);
% % % s4 = 1/(2*hx);
% % % 
% % % 
% % % Sigma_xx = zeros( n-2, n-2 );
% % % Sigma_xz = zeros( n-2, n-2 );
% % % Sigma_zz = zeros( n-2, n-2 );
% % % for p = 1:(n-2),
% % %     for q = 1:(n-2),
% % %         Sigma_xx( p, q ) = s1*( -u_sol( p, q+1 ) + u_sol( p+2, q+1 )) + s2*( -w_sol( p+1, q ) + w_sol( p+1, q+2 ));
% % %         Sigma_xz( p, q ) = s3*( -u_sol( p+1, q ) + u_sol( p+1, q+2 )) + s4*( -w_sol( p, q+1 ) + w_sol( p+2, q+1 ));
% % %         Sigma_zz( p, q ) = s1*( -w_sol( p, q+1 ) + w_sol( p+2, q+1 )) + s2*( -u_sol( p+1, q ) + u_sol( p+1, q+2 ));
% % %     end;
% % % end;
% % % 
% % % figure( 4 );
% % % surf( Sigma_xx );
% % % title( '\sigma_{xx} (Normal Stresses in X-direction)', 'FontSize', 36, 'FontWeight', 'bold');
% % % xlabel( 'X axis', 'FontSize', 30, 'FontWeight', 'bold' );
% % % ylabel( 'Z axis', 'FontSize', 30, 'FontWeight', 'bold' );
% % % zlabel( 'Normal Stress', 'FontSize', 30, 'FontWeight', 'bold' );
% % % 
% % % figure( 5 );
% % % surf( Sigma_zz );
% % % title( '\sigma_{zz} (Normal Stresses in Z-direction)', 'FontSize', 36, 'FontWeight', 'bold');
% % % xlabel( 'X axis', 'FontSize', 30, 'FontWeight', 'bold' );
% % % ylabel( 'Z axis', 'FontSize', 30, 'FontWeight', 'bold' );
% % % zlabel( 'Normal Stress', 'FontSize', 30, 'FontWeight', 'bold' );
% % % 
% % % figure( 6 );
% % % surf( Sigma_xz );
% % % title( '\sigma_{xz} = \sigma_{zx} (Shear Stresses in XZ-plane)', 'FontSize', 36, 'FontWeight', 'bold');
% % % xlabel( 'X axis', 'FontSize', 30, 'FontWeight', 'bold' );
% % % ylabel( 'Z axis', 'FontSize', 30, 'FontWeight', 'bold' );
% % % zlabel( 'Shear Stress', 'FontSize', 30, 'FontWeight', 'bold' );
toc;