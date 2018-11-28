function F=selfcons_anis(T,g,D,B,J0,M)

if nargin<6
    mode='selfcons';
end
if nargin<5
    plotdata='on';
end
if nargin<4
    Bmax=100;
    NB=100;
    B=[0 logspace(-3,log10(Bmax),NB)];
else
    NB=length(B);
    Bmax=max(B);
end
if nargin<3
    if strcmp(mode,'plot') 
        xdata=[];
    else
        xdata=0:0.1:10;
    end
end


B0=0:0.1:5;
[Mx2,Mz2]=anis_model_Spherical_2D_selfcons(T,[D,g,g],[],[0,0,0],B0,[],'on','selfcons');
%M2(isnan(M2))=0;
F=M2-M;

%{
figure(4321122)
cla
hold on
plot(B,M2,'k-')
plot(B,M,'r-')
drawnow
%}