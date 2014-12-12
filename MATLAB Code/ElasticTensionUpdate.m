function [dS_Stack] = ElasticTensionUpdate(R,S_Stack,sigma)

dS_Stack = zeros(3,1);

%   dS_Stack(1) = S_rr
%   dS_Stack(2) = z
%   dS_Stack(3) = r

if R==0
        dS_Stack(1)=0;
        dS_Stack(2)=(S_Stack(1)*d2corneal(0));
        dS_Stack(3)=(1/(S_Stack(1)*(1-sigma)+1));
else
        dS_Stack(1) = ( -(S_Stack(1)/R)*(1 + dcorneal(R)^2) ...
            + S_Stack(1)*dcorneal(R)*d2corneal(R)/(1 + dcorneal(R)^2) ...
            - dcorneal(R)*(sqrt(1 + dcorneal(R)^2)/R)*S_Stack(2) ...
            + (sigma*S_Stack(1) + tau(S_Stack(3))*(R - S_Stack(3))/S_Stack(3))*((1 + dcorneal(R)^2)/R));
        dS_Stack(2) = ( S_Stack(1)*dcorneal(R)/(R*sqrt(1 + dcorneal(R)^2)) );
        dS_Stack(3) = ( sqrt(1 + dcorneal(R)^2)*S_Stack(3)/( sqrt(1 + drest(S_Stack(3))^2)*( (1 + (1 - sigma^2)*S_Stack(1))*S_Stack(3) - sigma*((R-S_Stack(3))) ) ) );
end