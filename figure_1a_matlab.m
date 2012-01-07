clear all;

N = 5;

tic
M = importdata('figure_1a_0.1.dat',' ');

PP = M(:,1:N+1);
mu = M(:,N+2);
rho = M(:,N+3);

avnum = 100;

for i=1:length(M)/10
    mmu(i) = mean(mu((10*(i-1)+1):(10*i)));
    rrho(i) = mean(rho((10*(i-1)+1):(10*i)));
end

[B,IX] = sort(mmu);

plot(B,rrho(IX),'r');
hold on

M = importdata('figure_1a_0.3.dat',' ');

PP = M(:,1:N+1);
mu = M(:,N+2);
rho = M(:,N+3);

avnum = 100;

for i=1:length(M)/10
    mmu(i) = mean(mu((10*(i-1)+1):(10*i)));
    rrho(i) = mean(rho((10*(i-1)+1):(10*i)));
end

[B,IX] = sort(mmu);

plot(B,rrho(IX),'g');

M = importdata('figure_1a_0.5.dat',' ');

PP = M(:,1:N+1);
mu = M(:,N+2);
rho = M(:,N+3);

avnum = 100;

for i=1:length(M)/10
    mmu(i) = mean(mu((10*(i-1)+1):(10*i)));
    rrho(i) = mean(rho((10*(i-1)+1):(10*i)));
end

[B,IX] = sort(mmu);

plot(B,rrho(IX),'b');