function WVx = WVxToXY(xf,yf)
vx = 1;

if(abs(yf)<1e-12)
  w=0;
  t=xf/vx;
else
  w = 2*yf/(vx*(xf*xf+yf*yf));
  t = atan2(w*xf/vx, 1-w*yf/vx)/w;
end
    
 WVx = [w*t; vx*t];   
end