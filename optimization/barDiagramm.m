%x = [8 12 16 20 30 32 43 55 70 104];
x = [8; 12; 16; 20];
vals = [15.94 14.8 15.8 15.79; 14.35 13.88 13.89 14.46; 14.33 13.29 13.01 13.21; 13.54 13.0 12.63 13.44];
% y = [1079.8 629.82 534.16 338.15 249.38 184.11 151.38 103.84 79.51 76.51];
figure
% p = plot(x,y,'-r');
% title('Optimization progress as a function of active catheters');
% xlabel('Number of active catheters'); 
% ylabel('Fitness value'); 
b = bar(x,vals);
title('D_5 value of Rectum');
xlabel('Number of Channels'); 
ylabel('Dose Value');
legend('10 Generations','15 Generations', '20 Generations', 'Benchmark');