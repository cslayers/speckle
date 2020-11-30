
% drag the csv file into matlab comand line window
% then run the following code 
apple = sample;

pu = 0;
pv = 0;
pux = -0.11;
puy = -0.12;
pvx = -0.13;
pvy = -0.14;
puxx = 0.00011;
puxy = 0.00012;
puyy = 0.00013;
pvxx = 0.00014;
pvxy = 0.00015;
pvyy = 0.00016;

x = table2array(apple(:,2));
y = table2array(apple(:,3));

u = table2array(apple(:,4));
v = table2array(apple(:,7));
ux = table2array(apple(:,5));
uy = table2array(apple(:,6));
vx = table2array(apple(:,8));
vy = table2array(apple(:,9));
uxx = table2array(apple(:,25));
uxy = table2array(apple(:,26));
uyy = table2array(apple(:,27));
vxx = table2array(apple(:,28));
vxy = table2array(apple(:,29));
vyy = table2array(apple(:,30));

converge = table2array(apple(:,18));

disp(['ux mean:',num2str(mean(ux)),' std:',num2str(std(ux))]);
disp(['uy mean:',num2str(mean(uy)),' std:',num2str(std(uy))]);
disp(['vx mean:',num2str(mean(vx)),' std:',num2str(std(vx))]);
disp(['vy mean:',num2str(mean(vy)),' std:',num2str(std(vy))]);

disp(['uxx mean:',num2str(mean(uxx)),' std:',num2str(std(uxx))]);
disp(['uxy mean:',num2str(mean(uxy)),' std:',num2str(std(uxy))]);
disp(['uyy mean:',num2str(mean(uyy)),' std:',num2str(std(uyy))]);
disp(['vxx mean:',num2str(mean(vxx)),' std:',num2str(std(vxx))]);
disp(['vxy mean:',num2str(mean(vxy)),' std:',num2str(std(vxy))]);
disp(['vyy mean:',num2str(mean(vyy)),' std:',num2str(std(vyy))]);




x = x - 256;
y = y - 256;
xx = x .* x;
xy = x .* y;
yy = y .* y;

xp = x + pu + x * pux + y * puy + xx * puxx * 0.5 + xy * puxy + yy * puyy * 0.5;
yp = y + pv + x * pvx + y * pvy + xx * pvxx * 0.5 + xy * pvxy + yy * pvyy * 0.5;

erru = xp - x - u;
errv = yp - y - v;


[m n] = size(apple);

mbeu = sum(erru) / m;
mbev = sum(errv) / m;

sdeu = sqrt( sum(erru .* erru) / (m - 1) );
sdev = sqrt( sum(errv .* errv) / (m - 1) );

disp(['U MBE: ',num2str(mbeu)]);
disp(['V MBE: ',num2str(mbev)]);
disp(['U SDE: ',num2str(sdeu)]);
disp(['V SDE: ',num2str(sdev)]);


