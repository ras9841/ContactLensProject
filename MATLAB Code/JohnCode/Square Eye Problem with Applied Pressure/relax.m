dx=lx/nx;
dz=lz/nz;
x=dx*(0:nx);
z=dz*(0:nz);
%u=zeros(nx+1,nz+1);
%w=zeros(nx+1,nz+1);
%d=zeros(nx+1,nz+1);
sigma=.49;
k=1.0/(1-sigma);


for nit=1:1000
    for i=2:nx
        for j=2:nz
            d(i,j)=(u(i+1,j)-u(i-1,j))/(2*dx)+(w(i,j+1)-w(i,j-1))/(2*dz);
        end
    end
    for i=2:nx
        j=nz+1;
        w(i,j)=w(i,j-1)+dz*((1-sigma)*pressure(x(i)) - sigma*(u(i+1,j)-u(i-1,j))/(2*dx));
        u(i,j)=u(i,j-1)-dz*((w(i+1,j)-w(i-1,j))/(2*dx));
    end
    
    
        for i=2:nx

            for j=2:nz
               ures=(u(i+1,j)-2*u(i,j)+u(i-1,j))/(dx*dx) + (u(i,j+1)-2*u(i,j)+u(i,j-1))/(dz*dz) + k*(d(i+1,j)-d(i-1,j))/(2*dx);
               wres=(w(i+1,j)-2*w(i,j)+w(i-1,j))/(dx*dx) + (w(i,j+1)-2*w(i,j)+w(i,j-1))/(dz*dz) + k*(d(i,j+1)-d(i,j-1))/(2*dz);
               u(i,j)=u(i,j)+.2*dx*dx*ures;
               w(i,j)=w(i,j)+.2*dx*dx*wres;
            end
        end
    
        
        
end

