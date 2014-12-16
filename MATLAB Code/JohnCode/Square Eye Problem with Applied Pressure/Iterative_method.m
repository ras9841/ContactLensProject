x_iter = zeros(2*n*n,1);

for i=1:1000
   x_iter = x_iter + 0.001*(A_diag_slash_b - A_diag_slash_A*x_iter);   
   %norm(x_iter)
  % pause(1)
end

uw2 = reshape( x_iter, n, 2*n );
u_sol = uw2( 1:n, 1:n );
w_sol = uw2( 1:n, (n+1):2*n );

X = linspace(a,b,n);
Z = linspace(c,d,n);
[X,Z]=meshgrid(X,Z);

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
