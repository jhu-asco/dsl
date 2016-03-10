function f = car_ik()
% computes closed-form inverse kinematics solution for flows along
% SE(2) vector fields suitable for car-like trajectory generation
% Author: Marin Kobilarov

syms w1 v1 v2 w3 v3 t1 t2 t3 xf yf af real

syms k1 k2 real

pf = [xf;yf];

g2 = [1 0 t2*v2; 0 1 0; 0 0 1];

g = simplify(se2_exp(t1*[w1; v1; 0])*g2*se2_exp(t3*[w3; v3; 0]))

p = g(1:2,3);

p = simplify(subs(p, t3, (af - t1*w1)/w3 ))

e = p - pf

e = subs(e, {cos(t1*w1), sin(t1*w1)}, {(1-k1^2)/(1+k1^2), (2*k1)/(1+k1^2)})

[t2,k1] = solve(e, t2,k1)

simplify(t2)
simplify(k1)

% t1 = 2*atan(k1)/s1/w;
% t3 = (af/w - s1*t1)/s3;

% Run some test cases

w1=.5;
v1=1;
v2=1;
w3=.5;
v3=1;

for j=1:1

xf = rand*10-5; yf = rand*10-5; af = rand*2*pi-pi;

if j==1
  xf = 0; yf = 1; af = 0.1;
end

g = se2_g([af; xf; yf]);

% iterate through possible combos of signs of angular velocities
for k=1:3
  
  s2=[0;v2;0];
  switch k
   case 1,
    s1=[w1;v1;0];
    s3=[w3;v3;0];
   case 2,
    s1=[-w1;v1;0];
    s3=[w3;v3;0];
   case 3,
    s1=[w1;v1;0];
    s3=[-w3;v3;0];
   case 4,
    s1=[-w1;v1;0];
    s3=[-w3;v3;0];
   otherwise,
  end
  
[t1, t2, t3] = rs_times(s1, s2, s3, g)
if isinf(t1)
  continue
end

g=rs_evolve(s1, s2, s3, t1(1),t2(1),t3(1));
g=rs_evolve(s1, s2 ,s3, t1(2),t2(2),t3(2));

dt = .03;

g0 = eye(3);
se2_traj_plot(g0);

for i=1:2
g1s = se2_traj(g0, s1, t1(i), dt);
g2s = se2_traj(g1s(:,:,end), s2, t2(i), dt);
g3s = se2_traj(g2s(:,:,end), s3, t3(i), dt);
 
gs = cat(3, g0, g1s, g2s, g3s);
se2_traj_plot(gs);
axis equal
hold on
end


end
end

function [t1, t2,t3] = rs_times(s1,s2,s3,g)

w1 = s1(1);
v1 = s1(2);
w2 = s2(1);
v2 = s2(2);
w3 = s3(1);
v3 = s3(2);
xf = g(1,3); 
yf = g(2,3); 
ca = g(1,1); 
sa = g(2,1); 

q = (w1*w3*(2*v1*v3 - 2*v1*v3*ca + w1*w3*xf^2 + w1*w3*yf^2 - 2*v1*w3*yf ...
            + 2*v3*w1*yf*ca - 2*v3*w1*xf*sa));
if q < 0
  t1 = inf;
  t2 = inf;
  t3 = inf;
  return
end

t2 = [ q^(1/2)/(v2*w1*w3), -q^(1/2)/(v2*w1*w3)]; 
 
 
k1= [ ((w1*w3*(2*v1*v3 - 2*v1*v3*ca + w1*w3*xf^2 + w1*w3*yf^2 - 2*v1*w3*yf + 2*v3*w1*yf*ca - 2*v3*w1*xf*sa))^(1/2) + v3*w1*sa - w1*w3*xf)/(v3*w1 - 2*v1*w3 + v3*w1*ca + w1*w3*yf),...
 -((w1*w3*(2*v1*v3 - 2*v1*v3*ca + w1*w3*xf^2 + w1*w3*yf^2 - 2*v1*w3*yf + 2*v3*w1*yf*ca - 2*v3*w1*xf*sa))^(1/2) - v3*w1*sa + w1*w3*xf)/(v3*w1 - 2*v1*w3 + v3*w1*ca + w1*w3*yf)];
 
af = atan2(sa, ca);

t1 = 2*atan(k1)/w1;
af - t1*w1
w3
t3 = [fangle(af - t1(1)*w1), fangle(af - t1(2)*w1)]/w3


function g = rs_evolve(s1,s2,s3,t1,t2,t3)

w1 = s1(1);
v1 = s1(2);
w2 = s2(1);
v2 = s2(2);
w3 = s3(1);
v3 = s3(2);

g2 = [1 0 t2*v2; 0 1 0; 0 0 1];

g = se2_exp(t1*[w1; v1; 0])*g2*se2_exp(t3*[w3; v3; 0]);


function gs = se2_traj(g0, s, tf, dt)
if tf < 0
  dt = -dt;
end

N = floor(tf/dt);
gs = zeros(3,3,N);
t=0;
for i=1:N,
  t = i*dt;
  gs(:,:,i) = g0*se2_exp(t*s);
end
if (abs(t) < abs(tf))
  gs(:,:,N+1) = g0*se2_exp(tf*s);
end

function f = se2_traj_plot(gs)

ps = reshape(gs(1:2,3,:), 2,size(gs,3));

f = plot(ps(1,:), ps(2,:), '-');
hold on
quiver(ps(1,end), ps(2,end), gs(1,1,end), gs(2,1,end));

function a = fangle(a)
a = mod(a, 2*pi);
if a > pi
  a = a - 2*pi;
else
  if a < -pi
    a = a + 2*pi;
  end
end