R=load('R.csv');
W=load('W.csv');

figure
surf(R)
title ("R");

figure
surf(W)
title("W");

input('Press any key to close >> ');
