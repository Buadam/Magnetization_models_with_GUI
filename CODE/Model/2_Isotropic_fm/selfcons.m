function F=selfcons(T,g,S,B,J0,M)
muB=0.6720345504;


M2=Brillouin(T,[g,S],[],[0,0],B+2/g^2*J0/muB*M);
F=M2-M;

%{
figure(4321122)
cla
hold on
plot(B,M2,'k-')
plot(B,M,'r-')
drawnow
%}