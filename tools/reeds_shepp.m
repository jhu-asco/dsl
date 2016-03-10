function f = reeds_shepp()

syms v w t1 t2 t3 xf yf af s1 s3 real

syms k1 k2 real

pf = [xf;yf];

xi = [w; v; 0];

g2 = [1 0 t2*v; 0 1 0; 0 0 1];

g = simplify(se2_exp(t1*[s1*w; v; 0])*g2*se2_exp(t3*[s3*w; v; 0]))

p = g(1:2,3);

p = simplify(subs(p, t3, (af/w - s1*t1)/s3 ))

e = p - pf

e = subs(e, {cos(s1*t1*w), sin(s1*t1*w)}, {(1-k1^2)/(1+k1^2), (2*k1)/(1+k1^2)})

[t2,k1] = solve(e, t2,k1)

simplify(t2)
simplify(k1)

% t1 = 2*atan(k1)/s1/w;
% t3 = (af/w - s1*t1)/s3;



s1 = 1;
s3=1;
w=.5;
v=1;
[t1, t2, t3] = rs_times(v, w, s1, s3, 5, 5, pi/3)

g=rs_evolve(v,w,s1,s3,t1(1),t2(1),t3(1))

g=rs_evolve(v,w,s1,s3,t1(2),t2(2),t3(2))

for i=1:2

t = [t1(i), t2(i), t3(i)];
  

end 



function [t1, t2,t3] = rs_times(v,w,s1,s3,xf,yf,af)

t2 = [ (s1*s3*(2*v^2 - 2*v^2*cos(af) - 2*s3*v*w*yf + s1*s3*w^2*xf^2 + s1*s3*w^2*yf^2 - 2*s1*v*w*xf*sin(af) + 2*s1*v*w*yf*cos(af)))^(1/2)/(s1*s3*v*w), 
    -(s1*s3*(2*v^2 - 2*v^2*cos(af) - 2*s3*v*w*yf + s1*s3*w^2*xf^2 + s1*s3*w^2*yf^2 - 2*s1*v*w*xf*sin(af) + 2*s1*v*w*yf*cos(af)))^(1/2)/(s1*s3*v*w)];
 
k1= [ ((s1*s3*(2*v^2 - 2*v^2*cos(af) - 2*s3*v*w*yf + s1*s3*w^2*xf^2 + s1*s3*w^2*yf^2 - 2*s1*v*w*xf*sin(af) + 2*s1*v*w*yf*cos(af)))^(1/2) + s1*v*sin(af) - s1*s3*w*xf)/(s1*v - 2*s3*v + s1*v*cos(af) + s1*s3*w*yf),
 -((s1*s3*(2*v^2 - 2*v^2*cos(af) - 2*s3*v*w*yf + s1*s3*w^2*xf^2 + s1*s3*w^2*yf^2 - 2*s1*v*w*xf*             sin(af) + 2*s1*v*w*yf*cos(af)))^(1/2) - s1*v*sin(af) + s1*s3*w*xf)/(s1*v - 2*s3*v + s1*v*cos(af) + s1*s3*w*yf)];

t1 = 2*atan(k1)/s1/w;
t3 = (af/w - s1*t1)/s3;


function g = rs_evolve(v,w,s1,s3,t1,t2,t3)
xi = [w; v; 0];

g2 = [1 0 t2*v; 0 1 0; 0 0 1];

g = se2_exp(t1*[s1*w; v; 0])*g2*se2_exp(t3*[s3*w; v; 0]);
