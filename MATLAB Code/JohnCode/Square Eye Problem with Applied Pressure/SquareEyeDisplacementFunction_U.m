function f = SquareEyeDisplacementFunction_U( side, x, z )
switch side
    case 'left'
        f = 0;
    case 'right'
        f = 0; 
    case 'bottom'
        f = 0;
%     case 'top'
%         f = x;
    otherwise
        %y = 0;
        disp( 'Invalid displacement function argument of "side" in U.' )
end