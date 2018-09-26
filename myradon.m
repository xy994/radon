function [res,xp] = myradon(f,angle)

[N,M] = size(f);

xp=-92:1:92;
m = round(M/2);
n = round(N/2);

theta = angle:angle:180;
rhomax = ceil(sqrt(M^2+N^2));
rc = round(rhomax/2);
mt = length(theta);

res = cast(zeros(rhomax+3,mt),'double');

tic %?
for t = angle:angle:45
    costheta = cos(t*pi/180);
    sintheta = sin(t*pi/180);
    a = -costheta/sintheta;
    for r = 1:rhomax
        rho = r - rc;
        b = rho/sintheta;
        ymax = min(round(-a*m+b),n-1);
        ymin = max(round(a*m+b),-n);
        for y = ymin:ymax
            x = (y-b)/a;
            xfloor = floor(x);
            xup = x - xfloor;
            xlow = 1 - xup;
            x = xfloor;
            x = max(x, -m);
            x = min(x, m-2);
            res(rhomax-r+1, mt-find(theta==t)) = res(rhomax-r+1, mt-find(theta==t)) + xlow*f(y+n+1,x+m+1);
            res(rhomax-r+1, mt-find(theta==t)) = res(rhomax-r+1, mt-find(theta==t)) + xup*f(y+n+1,x+m+2);
        end
    end
end

for t = 45+angle:angle:90
    costheta = cos(t*pi/180);
    sintheta = sin(t*pi/180);
    a = -costheta/sintheta;
    for r = 1:rhomax
        rho = r - rc;
        b = rho/sintheta;
        xmax = min(round((-n-b)/a),m-1);
        xmin = max(round((n-b)/a),-m);
        for x = xmin:xmax
            y = a*x + b;
            yfloor = floor(y);
            yup = y - yfloor;
            ylow = 1 - yup;
            y = yfloor;
            y = max(y, -n);
            y = min(y, n-2);
            res(rhomax-r+1, mt-find(theta==t)) = res(rhomax-r+1, mt-find(theta==t)) + ylow*f(y+n+1,x+m+1);
            res(rhomax-r+1, mt-find(theta==t)) = res(rhomax-r+1, mt-find(theta==t)) + yup*f(y+n+2,x+m+1);
        end
    end
end

for t = 90+angle:angle:135
    costheta = cos(t*pi/180);
    sintheta = sin(t*pi/180);
    a = -costheta/sintheta;
    for r = 1:rhomax
        rho = r - rc;
        b = rho/sintheta;
        xmax = min(round((n-b)/a),m-1);
        xmin = max(round((-n-b)/a),-m);
        for x = xmin:xmax
            y = a*x + b;
            yfloor = floor(y);
            yup = y - yfloor;
            ylow = 1 - yup;
            y = yfloor;
            y = max(y, -n);
            y = min(y, n-2);
            res(rhomax-r+1, mt-find(theta==t)) = res(rhomax-r+1, mt-find(theta==t)) + ylow*f(y+n+1,x+m+1);
            res(rhomax-r+1, mt-find(theta==t)) = res(rhomax-r+1, mt-find(theta==t)) + yup*f(y+n+2,x+m+1);
        end
    end
end

for t = 135+angle:angle:180-angle
    costheta = cos(t*pi/180);
    sintheta = sin(t*pi/180);
    a = -costheta/sintheta;
    for r = 1:rhomax
        rho = r - rc;
        b = rho/sintheta;
        ymax = min(round(a*m+b),n-1);
        ymin = max(round(-a*m+b),-n);
        for y = ymin:ymax
            x = (y-b)/a;
            xfloor = floor(x);
            xup = x - xfloor;
            xlow = 1 - xup;
            x = xfloor;
            x = max(x, -m);
            x = min(x, m-2);
            res(rhomax-r+1, mt-find(theta==t)) = res(rhomax-r+1, mt-find(theta==t)) + xlow*f(y+n+1,x+m+1);
            res(rhomax-r+1, mt-find(theta==t)) = res(rhomax-r+1, mt-find(theta==t)) + xup*f(y+n+1,x+m+2);
        end
    end
end

for t = 180
    rhooffset = round((rhomax-M)/2);
    for x = 1:M
        r = x + rhooffset;
        r = rhomax - r +1;
        for y = 1:N
            res(r,find(theta==t)) = res(r,find(theta==t)) + f(y,x);
        end
    end
end

toc

rhoaxis = (1:rhomax+1)-rc;
figure
imagesc(1:180,rhoaxis,res);
colormap(hot), colorbar

end