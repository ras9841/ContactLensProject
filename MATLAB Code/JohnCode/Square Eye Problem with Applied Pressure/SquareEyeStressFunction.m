function f = SquareEyeStressFunction( x )
% x is the input argument 
% a,b  defines the interval over which x varies: x is on [a,b]
% for this stress function to be correct, and the displacements on the
% sides of the eye to be u(x,z) = w(x,z) = 0, then the stress function
% needs to have roots at a and b. If these conditions are not met, expect
% weird results at the top corners


f = -x*(1-x);%sin(2*pi*x);%-x*(1-x);

%  sigma = 0.1;
%  E = 1;
%  f = (-2*pi*E*sigma/((1+sigma)*(1-2*sigma)))*sin(2*pi*x);
end
