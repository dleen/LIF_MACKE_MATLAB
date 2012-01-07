clear all;

N = 5;

tic
M = importdata('figure_1a_0.1.dat',' ');

PP = M(:,1:N+1);
mu = M(:,N+2);
rho = M(:,N+3);

avnum = 5;

for i=1:length(M)/avnum
    mmu(i) = mean(mu((avnum*(i-1)+1):(avnum*i)));
    rrho(i) = mean(rho((avnum*(i-1)+1):(avnum*i)));
end

[B,IX] = sort(mmu);

plot(B,rrho(IX),'r');
hold on

M = importdata('figure_1a_0.3.dat',' ');

PP = M(:,1:N+1);
mu = M(:,N+2);
rho = M(:,N+3);

for i=1:length(M)/avnum
    mmu(i) = mean(mu((avnum*(i-1)+1):(avnum*i)));
    rrho(i) = mean(rho((avnum*(i-1)+1):(avnum*i)));
end

[B,IX] = sort(mmu);

plot(B,rrho(IX),'g');

M = importdata('figure_1a_0.5.dat',' ');

PP = M(:,1:N+1);
mu = M(:,N+2);
rho = M(:,N+3);

for i=1:length(M)/avnum
    mmu(i) = mean(mu((avnum*(i-1)+1):(avnum*i)));
    rrho(i) = mean(rho((avnum*(i-1)+1):(avnum*i)));
end

[B,IX] = sort(mmu);

plot(B,rrho(IX),'b');