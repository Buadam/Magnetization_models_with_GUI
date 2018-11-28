function Mav=Isotropic_AFM_mf(T,param_fit,param_fixed,fixed,xdata)

%% -Get parameters
if nargin<1
   T=2;
end
if nargin<2
    g=2; 
    S=5/2; 
    J0=-0.8;
else
    %Sort fixed and fitted parameters
    i=1;
    j=1;
    for l=1:length(fixed)
        if fixed(l)==1 %is the l-th parameter fixed?
            p(l)=param_fixed(i); 
            i=i+1;
        else
            p(l)=param_fit(j);
            j=j+1;
        end
    end
    g=p(1);
    S=p(2);
    J0=p(3);
end
if nargin<5
   xdata=0:0.1:10;
end

B1=xdata';
M_Br=Brillouin(T,[g,S],[],[0,0],B1);
M0=M_Br;

%MA=1.1*M0;
%MB=0.9*M0;

fun=@(M)selfcons_AFM(T,g,S,B1,J0,M);
options = optimoptions('fsolve');
M=fsolve(fun,[M0,-M0],options);
MA=M(:,1);
MB=M(:,2);
Mav=(MA+MB)/2;
%}

figure(312)
cla
hold all

plot(B1,MB,'g-.')
plot(B1,Mav,'k-')
plot(B1,MA,'b--')