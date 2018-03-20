plottig_surface(N)
load paramScan.dat
cl = zeros(N,N);
x = zeros(N,N);
y = zeros(N,N);
cl(:) = paramScan(:,3);
x(:) = paramScan(:,1);
y(:) = paramScan(:,2);
surf(x,y,cl)