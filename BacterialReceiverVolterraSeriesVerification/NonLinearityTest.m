% Non-Linearity Test


time = 0:0.01:10; 

x =@(a,t) a*t + a*(t.^2);

aArr = linspace(1,10,10);
figure
for i= 1:size(aArr,2)
subplot(3,4,i),plot(time, x(aArr(i),time), "r-o");
title("a:" + num2str(aArr(i)))
end