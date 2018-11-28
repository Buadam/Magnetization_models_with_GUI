function OUTPUT=fit_Anis2D_meanfield(T0,B0,param_fit,param_fixed,fixed,xdata,sizes,mode)
tic
%% -Get parameters
if nargin<8
    mode=[1,0,0,0];
end

if nargin<6
    xdata=linspace(0.01,Bmax,NB); 
end
if nargin<7
   sizes=[1,length(xdata)];
end
if nargin<5
    fixed=[0,0,0];
end
if nargin<4
    param_fixed=[];
end
if nargin<3
    gx=2; gz=2;
    D=10; J0=0.1; T0=2;
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
    D=p(1);
    gx=p(2);
    gz=p(3);
    J0=p(4);
end
if nargin<2
    B0=0.5;
end
if nargin<1
    T0=2;
end

%% - Define Polar grid
Thmax=pi/2;
NTh=61;
Th=linspace(0,Thmax,NTh);

OUTPUT=[];
if mode(1)==1
    [~,~,Mav]=Anis_2D_mf_mH(T0,[D,gx,gz,J0],xdata(sizes(1,1):sizes(1,2)),Th);
    OUTPUT=Mav;
end
if mode(2)==1
    [Mxx,~,~]=Anis_2D_mf_mH(T0,[D,gx,gz,J0],xdata(sizes(2,1):sizes(2,2)),Th);
    OUTPUT=[OUTPUT;Mxx];
end
if mode(3)==1
   [~,~,Mav_T]=Anis_2D_mf_mT(B0,[D,gx,gz,J0],xdata(sizes(3,1):sizes(3,2)),Th);
   OUTPUT=[OUTPUT;Mav_T];
end
if mode(4)==1
    [Mxx_T,~,~]=Anis_2D_mf_mT(B0,[D,gx,gz,J0],xdata(sizes(4,1):sizes(4,2)),Th);
    OUTPUT=[OUTPUT;Mxx_T];
end

toc

