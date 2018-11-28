function varargout = Anisotropic_paramagnet(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Anisotropic_paramagnet_OpeningFcn, ...
                   'gui_OutputFcn',  @Anisotropic_paramagnet_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%% - Initialize.
function Anisotropic_paramagnet_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

%Add Working directory to path
DataDir=fullfile(cd,'DATA');
WorkingDir = fullfile(cd,'CODE');
p = genpath(WorkingDir);
addpath(p)

%Clear graphs and variables
clear global;
cla(handles.mH_plot_2K)

%Read parameters
global D E gx gy gz S J0
D=str2double(get(handles.D,'String'));
E=str2double(get(handles.E,'String'));
gx=str2double(get(handles.gx,'String'));
gy=str2double(get(handles.gy,'String'));
gz=str2double(get(handles.gz,'String'));
S=str2double(get(handles.S,'String'));
J0=str2double(get(handles.J0,'String'));

Reset_fits(hObject, eventdata, handles);

%Load Measured Data
global mxH matlH mporH mporH5K mT mavT mxT 
mxH=dlmread(fullfile(DataDir,'mxH_2K_skal.dat'));   %tömeggel skálázott mx(H)@2K
matlH=dlmread(fullfile(DataDir,'matlH_2K_skal.dat'));     %tömeggel skálázott matl(H)@2K
mporH=dlmread(fullfile(DataDir,'mH_2K_por.dat'));    %porminta m(H)@2K (Sanya)

mporH5K=dlmread(fullfile(DataDir,'mH_5K_por.dat')); mporH5K(:,1)=mporH5K(:,1)/1e4;

mT=dlmread(fullfile(DataDir,'powder_mT.dat'));    %porminta m(T) 0.5T-ban
Bgr=dlmread(fullfile(DataDir,'CAP_M_T.txt'));     %háttér
mxT=dlmread(fullfile(DataDir,'mxT_skal.dat'));    %mxT tömeggel skálázva (normalas.m alapján)
mavT=dlmread(fullfile(DataDir,'matlT_skal.dat'));    %mxT tömeggel skálázva (normalas.m alapján)
mT(:,2)=mT(:,2)-interp1(Bgr(:,1),Bgr(:,2),mT(:,1),'linear','extrap'); %pormintából a háttér levonva

global mH mH2 mH3

mH=dlmread(fullfile(DataDir,'S1_Invivogen.dat'),'\t',0,0);
mH2=dlmread(fullfile(DataDir,'S2_SmallCrystals.dat'),'\t',0,0);
mH3=dlmread(fullfile(DataDir,'S3_BigCrystals.dat'),'\t',0,0);



%% -Plot.
function Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set(handles.mH_plot_2K,'nextplot','add')
%set(handles.mH_plot_2K,'xlabel','B(T)')
%set(handles.mH_plot_2K,'xlabel','M(\mu_B)')
global mxH matlH mporH mporH5K mT mavT mxT 
global mH mH2 mH3
global E D gx gy gz S J0

muB=9.274e-24; %[Am^2]
EMU=1e-3; %[Am^2]
M_hZ=4923.928/8;  %[g/mol]
conv=1e-3/muB*M_hZ/6e23; %from EMU to muB/fu.

%% - Define grids for (fast) plotting

%Magnetic field grid for m(H)
T=2;
Bmax=10;
NB=50;
%B=[0 logspace(-3,log10(Bmax),NB)];
B=linspace(0,Bmax,NB);
Thmax=pi/2;
NTh=61;
Theta=linspace(0,Thmax,NTh);

% Temperature-grid for m(T)

Tv=linspace(1,300,300);
B0=0.5;

switch get(handles.Model,'Value')
    %% -Isotropic Paramagnet
    case 1 
         disp('Model: Isotropic paramagnetic model (Brillouin(g,S))')
         
         % - plot mH
         disp('Model: Isotropic paramagnetic model (Brillouin(g,S))')
         Output = Brillouin(T,[gz,S],[],[0,0],B);
         % - plot mT
         Output_mT=Brillouin(Tv,[gz,S],[],[0,0],0.5);
         mx=Output; my=mx; mz=mx; mav=mx;  
         mxT_fit=Output_mT; 
         mzT_fit=mxT_fit;
         mavT_fit=mxT_fit;
         
         
    %% - Isotropic Heisenberg
    case 2
         disp('Model: Isotropic Heisenberg model with mean-field approximation')
         % - plot mH
         Output = Isotropic_mf_MH(T,[gz,S,J0],B);
         mx=Output; my=mx; mz=mx; mav=mx;  
         
         % - plot mT
         Output_mT=Isotropic_mf_MT(B0,[gz,S,J0],Tv);
         mxT_fit=Output_mT; 
         mzT_fit=mxT_fit;
         mavT_fit=mxT_fit;

         
    %% - Isotropic AFM Heisenberg     
    case 3 
         disp('Model: Isotropic AFM Heisenberg model with mean-field approximation')
         Output = Isotropic_AFM_mf(T,[gz,S,J0],[],[0,0,0],B);
         mx=Output; my=mx; mz=mx; mav=mx;      
    
    %% - Anisotropic Paramagnet 2D
    case 4 
        disp('Model: Paramagnet with axial anisotropy')
        [mx,mz,mav]=Anis_2D_mH(T,[D,gx,gz],B,Theta);
        [mxT_fit,mzT_fit,mavT_fit]=Anis_2D_mT(B0,[D,gx,gz],Tv,Theta);
    case 5 %Anisotropic Paramagnet 3D
        disp('Model: Paramagnet with triaxial anisotropy')
        Output = anis_model_Spherical_3D(T,[D,E,gx,gy,gz],[],[0,0,0,0,0],B,[],'off','plot');
        mx=Output(:,1); my=Output(:,2); mz=Output(:,3); mav=Output(:,4);
        mx5=Output(:,5); my5=Output(:,6); mz5=Output(:,7); mav5=Output(:,8);
        Tv=linspace(1,300,100);
        Output_mT = anis_model_Spherical_3D_mT(Tv,[D,E,gx,gy,gz],[],[0,0,0,0,0],B,[],'off','plot');
        mxT_fit=Output_mT(:,1); myT_fit=Output_mT(:,2); mzT_fit=Output_mT(:,3); mavT_fit=Output_mT(:,4);
        
    case 6 %Anisotropic Ferromagnet 2D
    disp('Model: Paramagnet with axial anisotropy and exchange interaction')
    [mx,mz,mav]=Anis_2D_mH(T,[D,gx,gz,J0],B,Theta);
    [mxT_fit,mzT_fit,mavT_fit]=Anis_2D_mT(B0,[D,gx,gz,J0],Tv,Theta);
    %{
    Output = Anisotropic_mf_single_pol(T,[D,gx,gz,J0],[],[0,0,0,0],B,[],[],'plot_Tdep');
    mx=Output{1};
    mz=Output{2};
    mav=Output{3};
    
        %mx=Output(:,1); my=mx; mz=Output(:,2); mav=Output(:,3);
        %mx5=Output(:,5); my5=Output(:,6); mz5=Output(:,7); mav5=Output(:,8);
        %Output_mT = anis_model_Spherical_3D_mT(Tv,[D,E,gx,gy,gz],[],[0,0,0,0,0],B,[],'off','plot');
     mxT_fit=Output{4}; 
     mzT_fit=Output{5};
     mavT_fit=Output{6};
    %}
end
    
%plot #1
plot(handles.mH_plot_2K,mxH(:,1)/1e4,mxH(:,2),'go')
set(handles.mH_plot_2K,'nextplot','add')
plot(handles.mH_plot_2K,matlH(:,1)/1e4,matlH(:,2),'bo')
plot(handles.mH_plot_2K,mporH(:,1),mporH(:,2),'kd')

%plot(handles.mH_plot_2K,mH(:,1)/1e4,mH(:,2)*conv,'kd')

plot(handles.mH_plot_2K,B,mx,'g-')
%plot(handles.mH_plot_2K,B,my,'g--')
plot(handles.mH_plot_2K,B,mz,'r-')
plot(handles.mH_plot_2K,B,mav,'b-')
set(handles.mH_plot_2K,'nextplot','replace')

%plot #3
plot(handles.mT_plot,mT(:,1),mT(:,2),'bo')
set(handles.mT_plot,'nextplot','add')
plot(handles.mT_plot,mxT(:,1),mxT(:,2),'go')

plot(handles.mT_plot,Tv,mxT_fit','g-')
plot(handles.mT_plot,Tv,mzT_fit','r-')
plot(handles.mT_plot,Tv,mavT_fit','b-')
set(handles.mT_plot,'nextplot','replace')

plot(handles.mT_inset,1./mT(:,1),mT(:,2),'bo')
set(handles.mT_inset,'nextplot','add')
plot(handles.mT_inset,1./mxT(:,1),mxT(:,2),'go')
plot(handles.mT_inset,1./Tv,mxT_fit','g-')
plot(handles.mT_inset,1./Tv,mzT_fit','r-')
plot(handles.mT_inset,1./Tv,mavT_fit','b-')
set(handles.mT_inset,'xlim',[0,0.25])

set(handles.mT_inset,'nextplot','replace')


%plot in separate windows
figure(1)
cla
hold all
box on
grid on
plot(mxH(:,1)/1e4,mxH(:,2),'go')
plot(matlH(:,1)/1e4,matlH(:,2),'bo')
plot(mporH(:,1),mporH(:,2),'kd')
plot(B,mx,'g-')
plot(B,mz,'r-')
plot(B,mav,'b-')
xlabel('B(T)')
ylabel('M(\mu_B/Fe^{3+})')

figure(3)
cla
hold all
box on
h1=plot(mT(:,1),mT(:,2),'bo');
plot(mxT(:,1),mxT(:,2),'go')
plot(Tv,mxT_fit','g-')
plot(Tv,mzT_fit','r-')
plot(Tv,mavT_fit','b-')
xlabel('T(K)')
ylabel('M(\mu_B/Fe^{3+})')


axes('Position',[.27 .32 .6 .6])
cla
hold on
box on
plot(1./mT(:,1),mT(:,2),'bo')
plot(1./mxT(:,1),mxT(:,2),'go')
plot(1./Tv,mxT_fit','g-')
plot(1./Tv,mzT_fit','r-')
plot(1./Tv,mavT_fit','b-')
xlim([0,0.25])
xlabel('1/T(K^{-1})')
ylabel('M(\mu_B/Fe^{3+})')


%{
%plot #2

plot(handles.mH_plot_5K,mporH5K(:,1),mporH5K(:,2),'kd')
set(handles.mH_plot_5K,'nextplot','add')

plot(handles.mH_plot_5K,B,mx5,'b-')
plot(handles.mH_plot_5K,B,my5,'g--')
plot(handles.mH_plot_5K,B,mz5,'r-')
plot(handles.mH_plot_5K,B,mav5,'k-')
set(handles.mH_plot_5K,'nextplot','replace')



%}
%Params for dataset#3: T=5K, D=0.25, E=0.2


%% - Fit selected model to data
function Fit_Callback(hObject, eventdata, handles)

global mxH matlH mporH mporH5K mT mavT mxT 
global mH mH2 mH3
global E D gx gy gz S J0

muB=9.274e-24; %[Am^2]
EMU=1e-3; %[Am^2]
M_hZ=4923.928/8;  %[g/mol]
conv=1e-3/muB*M_hZ/6e23; %from EMU to muB/fu.

xdata=[];
ydata=[];
mode=[0,0,0,0];
sizes=zeros(4,2);
if get(handles.fit_av,'Value')==1
    sizes(1,:)=[1,length(matlH(:,1))];
    xdata=[matlH(:,1)/1e4];
    ydata=[matlH(:,2)];
    mode(1)=1;
end
if get(handles.fit_ea,'Value')==1 
    sizes(2,:)=[length(xdata)+1;length(xdata)+length(mxH(:,1))];
    xdata=[xdata;mxH(:,1)/1e4];
    ydata=[ydata;mxH(:,2)];
    mode(2)=1;
end
if get(handles.fit_mT_av,'Value')==1
    sizes(3,:)=[length(xdata)+1;length(xdata)+length(mT(:,1))];
    xdata=[xdata;mT(:,1)];
    ydata=[ydata;mT(:,2)];
    mode(3)=1;
end
if get(handles.fit_mT_ea,'Value')==1
    sizes(4,:)=[length(xdata)+1;length(xdata)+length(mxT(:,1))];
    xdata=[xdata;mxT(:,1)];
    ydata=[ydata;mxT(:,2)];
    mode(4)=1;
end

Bmax=10;
NB=50;
B=[0 logspace(-3,log10(Bmax),NB)];
T=2;
B0=0.5;
switch get(handles.Model,'Value')
    case 1 %Isotropic Paramagnet
     
        fixed=[get(handles.Fix_gz,'Value'),get(handles.Fix_S,'Value')];
        param=[str2double(get(handles.gz,'String')),str2double(get(handles.S,'String'))];
        param0=[];
        param_fit=[];
        param_fix=[];
        j=1;
        k=1;
        for i=1:length(fixed)
            if fixed(i)==1
                param_fix(j)=param(i);
                j=j+1;
            else
                param_fit(k)=param(i);
                k=k+1;
                param0(i)=param(i);
            end
        end

        disp('Fit Model: Isotropic paramagnetic model (Brillouin(g,S))')
        Brfn=@(param_fit,xdata)fit_Brillouin(T,B0,param_fit,param_fix,fixed,xdata,sizes,mode); %Brilloin's function
        [param_fit,R,J,CovB] = nlinfit(xdata,ydata,Brfn,param0); %non-linear fitting yields the fitted parameters, the residuals, and the covariance matrix 
        [yFit, delta] = nlpredci(Brfn,xdata,param_fit,R,'cov',CovB); %predicting confidence intervals

        disp(['Fitted parameters: ' num2str(param_fit)])
        Confidence_var = nlparci(param_fit,R,'cov',CovB); %confidence intervals for the fitted parameters 
        Perror_var= sqrt(diag(CovB)); %error values estimated for the fitted parameters
        sprintf('Confidence Intervals: \n %f \t %f \n',Confidence_var')
        sprintf('Parameter Errors: \n %f \n %f \n  ',Perror_var')
        i=1;
        j=1;
        for l=1:length(fixed)
            if fixed(l)==1 %is the l-th parameter fixed?
                p(l)=param_fix(i); p_min(l)=p(l); p_max(l)=p(l);
                i=i+1;
            else
                p(l)=param_fit(j); p_min(l)=Confidence_var(j,1); p_max(l)=Confidence_var(j,2);
                j=j+1;
            end
        end

        Reset_fits(hObject, eventdata, handles);
        gz=p(1); gz_min=p_min(1); gz_max=p_max(1);
        S=p(2);  S_min=p_min(2); S_max=p_max(2);
        
        set(handles.gz_fit,'String',num2str(gz)); set(handles.gz_min,'String',num2str(gz_min)); set(handles.gz_max,'String',num2str(gz_max))
        set(handles.S_fit,'String',num2str(S)); set(handles.S_min,'String',num2str(S_min)); set(handles.S_max,'String',num2str(S_max))
        Plot_Callback(hObject, eventdata, handles)

    case 2 %Isotropic Heisenberg
        fixed=[get(handles.Fix_gz,'Value'),get(handles.Fix_S,'Value'),get(handles.Fix_J0,'Value')];
        param=[str2double(get(handles.gz,'String')),str2double(get(handles.S,'String')),str2double(get(handles.J0,'String'))];
        param0=[];
        param_fit=[];
        param_fix=[];
        j=1;
        k=1;
        for i=1:length(fixed)
            if fixed(i)==1
                param_fix(j)=param(i);
                j=j+1;
            else
                param_fit(k)=param(i);
                k=k+1;
                param0(i)=param(i);
            end
        end

        disp('Fit Model: Isotropic paramagnetic model (Brillouin(g,S))')
        modelfun=@(param_fit,xdata)fit_Isotropic_mf(T,B0,param_fit,param_fix,fixed,xdata,sizes,mode); %Brilloin's function

        [param_fit,R,J,CovB] = nlinfit(xdata,ydata,modelfun,param0); %non-linear fitting yields the fitted parameters, the residuals, and the covariance matrix 
        [yFit, delta] = nlpredci(modelfun,xdata,param_fit,R,'cov',CovB); %predicting confidence intervals

        disp(['Fitted parameters: ' num2str(param_fit)])
        Confidence_var = nlparci(param_fit,R,'cov',CovB); %confidence intervals for the fitted parameters 
        Perror_var= sqrt(diag(CovB)); %error values estimated for the fitted parameters
        sprintf('Confidence Intervals: \n %f \t %f \n',Confidence_var')
        sprintf('Parameter Errors: \n %f \n %f \n  ',Perror_var')
        i=1;
        j=1;
        for l=1:length(fixed)
            if fixed(l)==1 %is the l-th parameter fixed?
                p(l)=param_fix(i); p_min(l)=p(l); p_max(l)=p(l);
                i=i+1;
            else
                p(l)=param_fit(j); p_min(l)=Confidence_var(j,1); p_max(l)=Confidence_var(j,2);
                j=j+1;
            end
        end

        Reset_fits(hObject, eventdata, handles);
        gz=p(1); gz_min=p_min(1); gz_max=p_max(1);
        S=p(2);  S_min=p_min(2); S_max=p_max(2);
        J0=p(3);  J0_min=p_min(3); J0_max=p_max(3);
        
        set(handles.gz_fit,'String',num2str(gz)); set(handles.gz_min,'String',num2str(gz_min)); set(handles.gz_max,'String',num2str(gz_max))
        set(handles.S_fit,'String',num2str(S)); set(handles.S_min,'String',num2str(S_min)); set(handles.S_max,'String',num2str(S_max))
        set(handles.J0_fit,'String',num2str(J0)); set(handles.J0_min,'String',num2str(J0_min)); set(handles.J0_max,'String',num2str(J0_max))
        
        Plot_Callback(hObject, eventdata, handles)

    case 4 %Anisotropic Paramagnet 2D
        fixed=[get(handles.Fix_D,'Value'),get(handles.Fix_gx,'Value'),get(handles.Fix_gz,'Value')];
        param=[str2double(get(handles.D,'String')),str2double(get(handles.gx,'String')),str2double(get(handles.gz,'String'))];
        param0=[];
        param_fit=[];
        param_fix=[];
        j=1;
        k=1;
        for i=1:length(fixed)
            if fixed(i)==1
                param_fix(j)=param(i);
                j=j+1;
            else
                param_fit(k)=param(i);
                k=k+1;
                param0(i)=param(i);
            end
        end

        disp('Fit Model: Anisotropic paramagnetic model (D,gx,gz)')
        modelfun=@(param_fit,xdata)fit_anis2D(T,B0,param_fit,param_fix,fixed,xdata,sizes,mode); %Anisotropic paramagnet

        [param_fit,R,J,CovB] = nlinfit(xdata,ydata,modelfun,param0); %non-linear fitting yields the fitted parameters, the residuals, and the covariance matrix 
        [yFit, delta] = nlpredci(modelfun,xdata,param_fit,R,'cov',CovB); %predicting confidence intervals

        disp(['Fitted parameters: ' num2str(param_fit)])
        Confidence_var = nlparci(param_fit,R,'cov',CovB); %confidence intervals for the fitted parameters 
        Perror_var= sqrt(diag(CovB)); %error values estimated for the fitted parameters
        sprintf('Confidence Intervals: \n %f \t %f \n',Confidence_var')
        sprintf('Parameter Errors: \n %f \n %f \n  ',Perror_var')
        i=1;
        j=1;
        for l=1:length(fixed)
            if fixed(l)==1 %is the l-th parameter fixed?
                p(l)=param_fix(i); p_min(l)=p(l); p_max(l)=p(l);
                i=i+1;
            else
                p(l)=param_fit(j); p_min(l)=Confidence_var(j,1); p_max(l)=Confidence_var(j,2);
                j=j+1;
            end
        end

        Reset_fits(hObject, eventdata, handles);
        D=p(1); D_min=p_min(1); D_max=p_max(1);
        gx=p(2);  gx_min=p_min(2); gx_max=p_max(2);
        gz=p(3);  gz_min=p_min(3); gz_max=p_max(3);
        
        set(handles.D_fit,'String',num2str(D)); set(handles.D_min,'String',num2str(D_min)); set(handles.D_max,'String',num2str(D_max))
        set(handles.gx_fit,'String',num2str(gx)); set(handles.gx_min,'String',num2str(gx_min)); set(handles.gx_max,'String',num2str(gx_max))
        set(handles.gz_fit,'String',num2str(gz)); set(handles.gz_min,'String',num2str(gz_min)); set(handles.gz_max,'String',num2str(gz_max))
        
        Plot_Callback(hObject, eventdata, handles)


    case 5 %Anisotropic Paramagnet 3D
        fixed=[get(handles.Fix_D,'Value'),get(handles.Fix_E,'Value'),get(handles.Fix_gx,'Value'),get(handles.Fix_gy,'Value'),get(handles.Fix_gz,'Value')];
        param=[str2double(get(handles.D,'String')),str2double(get(handles.E,'String')),str2double(get(handles.gx,'String')),str2double(get(handles.gy,'String')),str2double(get(handles.gz,'String'))];
        param0=[];
        param_fit=[];
        param_fix=[];
        j=1;
        k=1;
        for i=1:length(fixed)
            if fixed(i)==1
                param_fix(j)=param(i);
                j=j+1;
            else
                param_fit(k)=param(i);
                k=k+1;
                param0(i)=param(i);
            end
        end

        disp('Fit Model: Anisotropic paramagnetic model (D,gx,gz)')
        modelfun=@(param_fit,xdata)anis_model_Spherical_3D(T,param_fit,param_fix,fixed,B,xdata,'off','fit'); %Brilloin's function

        [param_fit,R,J,CovB] = nlinfit(xdata,ydata,modelfun,param0); %non-linear fitting yields the fitted parameters, the residuals, and the covariance matrix 
        [yFit, delta] = nlpredci(modelfun,xdata,param_fit,R,'cov',CovB); %predicting confidence intervals

        disp(['Fitted parameters: ' num2str(param_fit)])
        Confidence_var = nlparci(param_fit,R,'cov',CovB); %confidence intervals for the fitted parameters 
        Perror_var= sqrt(diag(CovB)); %error values estimated for the fitted parameters
        sprintf('Confidence Intervals: \n %f \t %f \n',Confidence_var')
        sprintf('Parameter Errors: \n %f \n %f \n  ',Perror_var')
        i=1;
        j=1;
        for l=1:length(fixed)
            if fixed(l)==1 %is the l-th parameter fixed?
                p(l)=param_fix(i); p_min(l)=p(l); p_max(l)=p(l);
                i=i+1;
            else
                p(l)=param_fit(j); p_min(l)=Confidence_var(j,1); p_max(l)=Confidence_var(j,2);
                j=j+1;
            end
        end

        Reset_fits(hObject, eventdata, handles);
        D=p(1); D_min=p_min(1); D_max=p_max(1);
        E=p(2); E_min=p_min(2); E_max=p_max(2);
        gx=p(3);  gx_min=p_min(3); gx_max=p_max(3);
        gy=p(4);  gx_min=p_min(4); gx_max=p_max(4);
        gz=p(5);  gz_min=p_min(5); gz_max=p_max(5);
        
        set(handles.D_fit,'String',num2str(D)); set(handles.D_min,'String',num2str(D_min)); set(handles.D_max,'String',num2str(D_max))
        set(handles.E_fit,'String',num2str(E)); set(handles.E_min,'String',num2str(E_min)); set(handles.E_max,'String',num2str(E_max))
        set(handles.gx_fit,'String',num2str(gx)); set(handles.gx_min,'String',num2str(gx_min)); set(handles.gx_max,'String',num2str(gx_max))
        set(handles.gy_fit,'String',num2str(gy)); set(handles.gy_min,'String',num2str(gy_min)); set(handles.gy_max,'String',num2str(gy_max))
        set(handles.gz_fit,'String',num2str(gz)); set(handles.gz_min,'String',num2str(gz_min)); set(handles.gz_max,'String',num2str(gz_max))
        
        Plot_Callback(hObject, eventdata, handles)
        
        
        case 6 %Anisotropic Paramagnet 2D
        fixed=[get(handles.Fix_D,'Value'),get(handles.Fix_gx,'Value'),get(handles.Fix_gz,'Value'),get(handles.Fix_J0,'Value')];
        param=[str2double(get(handles.D,'String')),str2double(get(handles.gx,'String')),str2double(get(handles.gz,'String')),str2double(get(handles.J0,'String'))];
        param0=[];
        param_fit=[];
        param_fix=[];
        j=1;
        k=1;
        for i=1:length(fixed)
            if fixed(i)==1
                param_fix(j)=param(i);
                j=j+1;
            else
                param_fit(k)=param(i);
                k=k+1;
                param0(i)=param(i);
            end
        end

        disp('Fit Model: Anisotropic paramagnetic model with exchange (D,gx,gz,J0)')
        modelfun=@(param_fit,xdata)Anisotropic_mf_single_pol(T,param_fit,param_fix,fixed,B,xdata,sizes,'fit_av_ea_MT'); %Brilloin's function

        [param_fit,R,J,CovB] = nlinfit(xdata,ydata,modelfun,param0); %non-linear fitting yields the fitted parameters, the residuals, and the covariance matrix 
        [yFit, delta] = nlpredci(modelfun,xdata,param_fit,R,'cov',CovB); %predicting confidence intervals

        disp(['Fitted parameters: ' num2str(param_fit)])
        Confidence_var = nlparci(param_fit,R,'cov',CovB); %confidence intervals for the fitted parameters 
        Perror_var= sqrt(diag(CovB)); %error values estimated for the fitted parameters
        sprintf('Confidence Intervals: \n %f \t %f \n',Confidence_var')
        sprintf('Parameter Errors: \n %f \n %f \n  ',Perror_var')
        i=1;
        j=1;
        for l=1:length(fixed)
            if fixed(l)==1 %is the l-th parameter fixed?
                p(l)=param_fix(i); p_min(l)=p(l); p_max(l)=p(l);
                i=i+1;
            else
                p(l)=param_fit(j); p_min(l)=Confidence_var(j,1); p_max(l)=Confidence_var(j,2);
                j=j+1;
            end
        end

        Reset_fits(hObject, eventdata, handles);
        D=p(1); D_min=p_min(1); D_max=p_max(1);
        gx=p(2);  gx_min=p_min(2); gx_max=p_max(2);
        gz=p(3);  gz_min=p_min(3); gz_max=p_max(3);
        J0=p(4);  J0_min=p_min(4); J0_max=p_max(4);
        
        set(handles.D_fit,'String',num2str(D)); set(handles.D_min,'String',num2str(D_min)); set(handles.D_max,'String',num2str(D_max))
        set(handles.gx_fit,'String',num2str(gx)); set(handles.gx_min,'String',num2str(gx_min)); set(handles.gx_max,'String',num2str(gx_max))
        set(handles.gz_fit,'String',num2str(gz)); set(handles.gz_min,'String',num2str(gz_min)); set(handles.gz_max,'String',num2str(gz_max))
        set(handles.J0_fit,'String',num2str(J0)); set(handles.J0_min,'String',num2str(J0_min)); set(handles.J0_max,'String',num2str(J0_max))
        
        Plot_Callback(hObject, eventdata, handles)

        
        
end


% --- Outputs from this function are returned to the command line.
function varargout = Anisotropic_paramagnet_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%% --- Object callbacks
function D_Callback(hObject, eventdata, handles)
global D
D=str2double(get(handles.D,'string'));
set(handles.D_Slider,'Value',D);

if get(handles.Instantplot,'Value')==1
   Plot_Callback(hObject, eventdata, handles)
end
function E_Callback(hObject, eventdata, handles)
global E 
E=str2double(get(handles.E,'string'));
set(handles.E_Slider,'Value',E);
if get(handles.Instantplot,'Value')==1 
    Plot_Callback(hObject, eventdata, handles)
end

function gx_Callback(hObject, eventdata, handles)
global gx gy
gx=str2double(get(handles.gx,'string'));
    if get(handles.gx_gy,'Value')==1
        gy=gx;
        set(handles.gy_Slider,'Value',gy);
    end
set(handles.gx_Slider,'Value',gx);
if get(handles.Instantplot,'Value')==1   
    Plot_Callback(hObject, eventdata, handles) 
end

function gy_Callback(hObject, eventdata, handles)
global gx gy
gy=str2double(get(handles.gy,'string'));
    if get(handles.gx_gy,'Value')==1
        gx=gy;
        set(handles.gx_Slider,'Value',gx);
    end
set(handles.gy_Slider,'Value',gy);
if get(handles.Instantplot,'Value')==1 
    Plot_Callback(hObject, eventdata, handles)
end

function gz_Callback(hObject, eventdata, handles)
global gz
gz=str2double(get(handles.gz,'string'));
set(handles.gz_Slider,'Value',gz);
if get(handles.Instantplot,'Value')==1 
    Plot_Callback(hObject, eventdata, handles) 
end

function S_Callback(hObject, eventdata, handles)
global S
S=str2double(get(handles.S,'string'));
set(handles.S_Slider,'Value',S);
if get(handles.Instantplot,'Value')==1
    Plot_Callback(hObject, eventdata, handles)
end

function J0_Callback(hObject, eventdata, handles)
global J0
J0=str2double(get(handles.J0,'string'));
set(handles.J0_Slider,'Value',J0);

if get(handles.Instantplot,'Value')==1  
    Plot_Callback(hObject, eventdata, handles)
end

%% --- Executes on slider movement.
function D_Slider_Callback(hObject, eventdata, handles)
global D
D=get(handles.D_Slider,'value');
set(handles.D,'string',num2str(D));
if get(handles.Instantplot,'Value')==1  
    Plot_Callback(hObject, eventdata, handles)
end

function E_Slider_Callback(hObject, eventdata, handles)
global E
E=get(handles.E_Slider,'value');
set(handles.E,'string',num2str(E));
if get(handles.Instantplot,'Value')==1 
    Plot_Callback(hObject, eventdata, handles)
end

function gx_Slider_Callback(hObject, eventdata, handles)
global gx gy
gx=get(handles.gx_Slider,'value');
set(handles.gx,'string',num2str(gx));
if get(handles.gx_gy,'Value')==1
        gy=gx;
        set(handles.gy,'string',num2str(gy));
end
if get(handles.Instantplot,'Value')==1   
    Plot_Callback(hObject, eventdata, handles)
end

function gy_Slider_Callback(hObject, eventdata, handles)
global gx gy
gy=get(handles.gy_Slider,'value');
set(handles.gy,'string',num2str(gy));
if get(handles.gx_gy,'Value')==1
        gx=gy;
        set(handles.gx,'string',num2str(gx));
end
if get(handles.Instantplot,'Value')==1  
    Plot_Callback(hObject, eventdata, handles)
end

function gz_Slider_Callback(hObject, eventdata, handles)
global gz
gz=get(handles.gz_Slider,'value');
set(handles.gz,'string',num2str(gz));
if get(handles.Instantplot,'Value')==1 
    Plot_Callback(hObject, eventdata, handles)
end


function S_Slider_Callback(hObject, eventdata, handles)
global S
S=get(handles.S_Slider,'value');
set(handles.S,'string',num2str(S));
if get(handles.Instantplot,'Value')==1   
    Plot_Callback(hObject, eventdata, handles)
end

function J0_Slider_Callback(hObject, eventdata, handles)
global J0
J0=get(handles.J0_Slider,'value');
set(handles.J0,'string',num2str(J0));
if get(handles.Instantplot,'Value')==1 
    Plot_Callback(hObject, eventdata, handles)
end

%% --- Executes on selection change in Model.
function Model_Callback(hObject, eventdata, handles)
switch get(handles.Model,'Value')
    case 1 %Isotropic Paramagnet
            set(handles.gx_gy,'Value',1)
    case 2 %Isotropic Heisenberg
            set(handles.gx_gy,'Value',1)
    case 3 %Anisotropic Paramagnet 2D
            set(handles.gx_gy,'Value',1)
    case 4 %Anisotropic Paramagnet 3D
            set(handles.gx_gy,'Value',0)
end

% --- Executes on button press in LogX.
function LogX_Callback(hObject, eventdata, handles)
if get(handles.LogX,'value')==1
    set(handles.mH_plot_2K,'xscale','log')
else
    set(handles.mH_plot_2K,'xscale','lin')
end

% --- Executes on button press in LogY.
function LogY_Callback(hObject, eventdata, handles)
if get(handles.LogY,'value')==1
    set(handles.mH_plot_2K,'yscale','log')
else
    set(handles.mH_plot_2K,'yscale','lin')
end

% --- Executes on button press in Fix_gy.
function Fix_gy_Callback(hObject, eventdata, handles)
function Fix_gx_Callback(hObject, eventdata, handles)
function Fix_D_Callback(hObject, eventdata, handles)
function Fix_E_Callback(hObject, eventdata, handles)
function Fix_gz_Callback(hObject, eventdata, handles)
function Fix_J0_Callback(hObject, eventdata, handles)
function gx_gy_Callback(hObject, eventdata, handles)
function Fix_S_Callback(hObject, eventdata, handles)

%%
%Create functions
% --- Executes during object creation, after setting all properties.

function D_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function E_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gx_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gy_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gz_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function J0_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function D_Slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function E_Slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function gx_Slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function gy_Slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function gz_Slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function J0_Slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Model_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit7_Callback(hObject, eventdata, handles)


function edit7_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function S_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function S_Slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function D_fit_Callback(hObject, eventdata, handles)

function D_fit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function E_fit_Callback(hObject, eventdata, handles)

function E_fit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gx_fit_Callback(hObject, eventdata, handles)

function gx_fit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gy_fit_Callback(hObject, eventdata, handles)

function gy_fit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gz_fit_Callback(hObject, eventdata, handles)

function gz_fit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function J0_fit_CreateFcn(hObject, eventdata, handles)

function S_fit_CreateFcn(hObject, eventdata, handles)



function D_max_Callback(hObject, eventdata, handles)

function D_max_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function E_max_Callback(hObject, eventdata, handles)

function E_max_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gx_max_Callback(hObject, eventdata, handles)

function gx_max_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gy_max_Callback(hObject, eventdata, handles)

function gy_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gz_max_Callback(hObject, eventdata, handles)

function gz_max_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function D_min_Callback(hObject, eventdata, handles)

function D_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function E_min_Callback(hObject, eventdata, handles)

function E_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gx_min_Callback(hObject, eventdata, handles)
function gx_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gy_min_Callback(hObject, eventdata, handles)
function gy_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gz_min_Callback(hObject, eventdata, handles)

function gz_min_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Instantplot.
function Instantplot_Callback(hObject, eventdata, handles)
% hObject    handle to Instantplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Instantplot


% --- Executes on button press in fit_av.
function fit_av_Callback(hObject, eventdata, handles)
% hObject    handle to fit_av (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_av


% --- Executes on button press in fit_ea.
function fit_ea_Callback(hObject, eventdata, handles)
% hObject    handle to fit_ea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_ea


% --- Executes on button press in fit_mT_av.
function fit_mT_av_Callback(hObject, eventdata, handles)
% hObject    handle to fit_mT_av (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_mT_av


% --- Executes on button press in fit_mT_ea.
function fit_mT_ea_Callback(hObject, eventdata, handles)
% hObject    handle to fit_mT_ea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_mT_ea
