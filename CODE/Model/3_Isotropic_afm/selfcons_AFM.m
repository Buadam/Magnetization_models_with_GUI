function F=selfcons_AFM(T,g,S,B,J0,M)
muB=0.6720345504;

MA=M(:,1);
MB=M(:,2);
MA2=Brillouin(T,[g,S],[],[0,0],B+2/g^2*J0/muB*MB);
MB2=Brillouin(T,[g,S],[],[0,0],B+2/g^2*J0/muB*MA);
Br=Brillouin(T,[g,S],[],[0,0],B);

M2=[MA2,MB2];


figure(311)
cla
hold all
plot(B,MA,'r-')
plot(B,MA2,'b-')
plot(B,MA,'g-')
plot(B,MA2,'k-')
drawnow

figure(312)
cla
hold all
plot(B,(MA2+MB2)/2,'k-')
plot(B,MB2,'r-.')
plot(B,MA2,'b--')
plot(B,Br,'g-')
drawnow


F=M2-M;

