function dxdt = Lotka_Volterra(x, t)
dx1dt = 1.1*x(1) - 0.4*x(1)*x(2); 
dx2dt = 0.1*x(2)*x(1) - 0.4*x(2); 
dxdt = [dx1dt; dx2dt]; 
end