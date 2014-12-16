% Poisson's Ratio (must be on interval of [-1,1/2))
sigma = 0.25; % sigma = 0.4998 works well

n = 101; % number of POINTS in x-direction and z-direction for both 
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

% construct the computational grid
X = linspace(a,b,n);
Z = linspace(c,d,n);
[X,Z]=meshgrid(X,Z);


% construct pressure and the derivative of the pressure vector
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


% parameters for bcs
if ( k == k1 ) && ( k1 ~= k2 )
    alpha4 = ( E*sigma*hz ) / (2*(1-sigma^2)*(2-sigma)*hx^2 );
    alpha5 = ( -E ) / ( (1-sigma^2)*hz );
    alpha6 = ( E*sigma ) / ( 2*(1-sigma^2)*hx );
elseif k == k2
    alpha4 = ( E*sigma*hz ) / ( 2*(1+sigma)*(1-2*sigma)*hx^2 );
    alpha5 = ( -E*(1-sigma) ) / ( (1+sigma)*(1-2*sigma)*hz );
    alpha6 = ( E*sigma ) / ( 2*(1+sigma)*(1-2*sigma)*hx );
end;

if ( k == k1 ) && ( k1 ~= k2 )
    alpha1 = ( -2*E ) / ( (1+sigma)*hx^2 );
    alpha2 = ( -2*E ) / ( (1+sigma)*hz^2 );
    alpha3 = ( E ) / ( (1+sigma)*hx*hz );
elseif k == k2
    alpha1 = ( -E*(2-sigma) ) / ( (1+sigma)*hx^2 );
    alpha2 = ( -2*E*(1-sigma) ) / ( (1+sigma)*hz^2 );
    alpha3 = ( E*(1-sigma) ) / ( (1+sigma)*hx*hz );
end;


% construct initial u and w displacements

% u = zeros(n,n);
% w = zeros(n,n);
% u_new = u;
% w_new = w;
omega = .001;

% relaxation iterative loop

    for r=1:100 %relaxation loop
        for i=2:n %we exclude z=0 because u=w=0
            for j=2:n-1 %we exclude x=0 and x=1 because u=w=0 there
    
               if i==n
                   %boundary condition
                   u_new(i,j) = u(i,j) + omega*hx*hx*( alpha1*(u(i,j-1) -2*u(i,j) + u(i,j+1)) ...
                       + alpha2*(u(i-1,j)-u(i,j)) + alpha3*(w(i,j+1)-w(i,j-1)) - vec_p_prime(j));
                   w_new(i,j) = w(i,j) + omega*hz*hz*( alpha4*(w(i,j+1) -2*w(i,j)+ w(i,j-1)) ...
                       + alpha5*(w(i-1,j)-w(i,j))+ alpha6*(u(i,j+1)-u(i,j-1)) - vec_p(j));
               else
                  %inside the domain
                  u_new(i,j) = u(i,j) + 0.01*hx*hx*( ...
                      ((k+1)/hx^2)*(u(i,j+1)-2*u(i,j)+u(i,j-1)) + (1/hz^2)*(u(i+1,j)-2*u(i,j)+u(i-1,j)) ...
                      + (k/(4*hx*hz))*( w(i+1,j+1) -w(i+1,j-1) -w(i-1,j+1) +w(i-1,j-1) ));
                  
                  w_new(i,j) = w(i,j) + 0.01*hz*hz*( ...
                      (1/hx^2)*(w(i,j+1)-2*w(i,j)+w(i,j-1)) + ((k+1)/hz^2)*(w(i+1,j)-2*w(i,j)+w(i-1,j)) ...
                      + (k/(4*hx*hz))*( u(i+1,j+1) -u(i+1,j-1) -u(i-1,j+1) +u(i-1,j-1) ));
               end
            end
            u=u_new;
            w=w_new;
        end
    end
