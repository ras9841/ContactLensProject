function [value,isterminal,direction]=ElasticTensionShootingEvent(R,S_Stack,E,sigma,tau)

value(1) = S_Stack(3)-1;
isterminal(1) = 1;
direction(1) = 1;