load paramScan.dat
cl = zeros(1000,1000);
x = zeros(1000,1000);
y = zeros(1000,1000);
cl(:) = paramScan(:,3);
x(:) = paramScan(:,1);
y(:) = paramScan(:,2);
surf(x,y,cl)