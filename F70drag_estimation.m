clear all
close all
clc

%rohith8272@gmail.com


%%  drag estimation of Fokker 70
% This script computes the parasitic drag using component drag build up
% method and also estimation of the induced drag and the thrust line as a
% variation of velocity in mach number

%% environmental variables
v	    =1.48*10.^-5  ;   %kinematic viscositydensity by rho
rho     =0.52 ;           %density at 8000 m alt above MSL
rho0    =1.225  ;             %density at MSL
vcruise =200    ;              %cruise velocity 1
Vcruise2=220     ;           %cruise veloicty 2
u       =1.789*10.^-5  ;      %dynamic viscosity
a       =340.2778      ;        %airspeed at 8000 m
Vapp    =61.2189       ;         %F 70 approach speed in m/s
M       =vcruise/a;
%% parasitic drag component estimation

%% wetted area calculation

%gearpod
GP_length=5.5;
GP_width=3.3;
GP_height=0.75;
TC_GP=GP_height/GP_length;
Sw_GP=0;
Airfoil_area=(0.5*0.539*1.34)+(0.5*0.539*4.04);


GP_area_bottom=(GP_length*GP_width);
GP_area_front=(GP_height*GP_width)-(0.6286);
GP_area_side=(GP_length*GP_height)-(2*0.9*Airfoil_area)   ; %% 90% of wing root airfoil area shadows gearpod's side area

Swet_GP=GP_area_bottom+GP_area_front+GP_area_side;

% fuselage
L1=3.63 ;          %cockpit section
L2=14.95  ;        %passenger section
L3=9.2    ;        %rear upsweep section
Dfuse=3.3;
d1=3.3;
d2=0.245    ;       %averaged diametere at fuselage tip
fuse_length=27.88;
Fuse_width=3.3;

 Swet_fuse_1  =(((pi*Dfuse)*(12*L1.^2).^-1)*((4*L1.^2 + ((Dfuse.^2)/4)).^1.5-((Dfuse.^3)*8.^-1))-((pi*Dfuse.^2)/4));% praboloid  aprroximation
 Swet_fuse_2  =(pi*Dfuse*L2);%cylinder approximation
 Swet_fuse_3  =((pi*(d1+d2)*0.5)*sqrt((L3.^2)+(((d1.^2)-(d2.^2))/4)));  % frustum approximation
 Swet_fuse    = Swet_fuse_1+Swet_fuse_2+Swet_fuse_3;  % total fuselage wetted area

%wing
Sref     =93.5  ; %reference area for all calculations is the wing area
MAC_wing =3.8;
tc_wing  =0.102750;
Sw_wing= 17.5;
lambda_wing=0.235;

Swet_wing=(2+(0.1057/3))*Sref-(2*GP_area_bottom);% the airfoil is NACA 0010 with 10.57 pc thickness( ref anatomy of an airplane(Darrol stinton)

%horizontal tail
Sht=21.72;
lambda_ht=0.39;
Cr_ht=1.214;
MAC_ht=Cr_ht*(2/3)*((1+lambda_ht+(lambda_ht.^2)).*(1+lambda_ht).^-1);
Swet_ht=(2+(0.10/3))*Sht ; % for NACA 0010
Sw_ht=26;


%vertical tail
Svt=12.3;
lambda_vt=0.74;
Cr_vt=4.2842;
MAC_vt=Cr_vt*(2/3)*((1+lambda_vt+(lambda_vt.^2)).*(1+lambda_vt).^-1);
Swet_vt=(2+(0.10/3))*Svt;  % for NACA 0010
Sw_vt=41;
TC_HT=0.10;
TC_VT=0.10;

%nacelle
nacelle_length =5.13;
nacelle_length_half=nacelle_length/2;
nacelle_width  =1.7;
Swet_nacelle=87;
nacelle_D1=1.7;     %%Maximum nacelle diameter
nacelle_D2=0.9;     %% Nacelle frontal diameter

Swet_nacelle_frust=2*((pi*(nacelle_D1 + nacelle_D2)/2)*sqrt((nacelle_length^2/4)+(nacelle_D1^2-nacelle_D2^2)/4));  
Swet_nacelle=4*(Swet_nacelle_frust);



% pylon exposed dimensions
TC_pylon=0.15;
Cr_pylon=4.74;
Ct_pylon=4.0;
Le=0.41;
Te=1.35;
ar_py1=Cr_pylon*Le;
tr_base=0.94;
ar_py2=0.5*Cr_pylon*tr_base;
S_pylon=ar_py1+ar_py2;
Swet_pylon=2*((2+(TC_pylon/3))*S_pylon);

MAC_pylon=4.36406;
Sw_pylon=0;



Swet=[Swet_wing,Swet_fuse,Swet_ht,Swet_vt,Swet_pylon,Swet_nacelle] ;                                       

%reynolds number calculation

Re_wing=(rho*vcruise*MAC_wing)*u.^-1;
Re_HT=(rho*vcruise*MAC_ht)*u.^-1;
Re_VT=(rho*vcruise*MAC_vt)*u.^-1;
Re_pylon=(rho*vcruise*MAC_pylon)*u.^-1;
Re_fuselage=(rho*vcruise*fuse_length)*u.^-1;
Re_nacelle=(rho*vcruise*nacelle_length)*u.^-1;
Re_gearpod=(rho*vcruise*GP_length)*u.^-1;
Re=[Re_wing;Re_fuselage;Re_HT;Re_VT;Re_pylon;Re_nacelle;Re_gearpod];


V=Vapp:1:340
machsweep=V/a

%skin friciton calculation

 %cf(zeros5,1)
Cf_wing    =0.455*((log10(Re_wing).^2.58)*(1+(0.144*machsweep.^2))).^-1;
Cf_fuselage=0.455*((log10(Re_fuselage).^2.58)*(1+(0.144*machsweep.^2))).^-1;
Cf_HT      =0.455*((log10(Re_HT).^2.58)*(1+(0.144*machsweep.^2))).^-1;
Cf_VT      =0.455*((log10(Re_VT).^2.58)*(1+(0.144*machsweep.^2))).^-1;
Cf_nacelle =0.455*((log10(Re_nacelle).^2.58)*(1+(0.144*machsweep.^2))).^-1;
Cf_pylon   =0.455*((log10(Re_pylon).^2.58)*(1+(0.144*machsweep.^2))).^-1;
Cf_gearpod =0.455*((log10(Re_gearpod).^2.58)*(1+(0.144*machsweep.^2))).^-1;

%Cf_turb    =0.455*((log(Re).^2.58)*(1+(0.144*machsweep.^2))).^-1
Cf         =[Cf_wing;Cf_fuselage;Cf_HT;Cf_VT;Cf_pylon;Cf_nacelle;Cf_gearpod];



%form factor calculation
xc_wing=0.25*MAC_wing;
xc_Htail=0.25*MAC_ht;
xc_Vtail=0.25*MAC_vt;
xc_pylon=0.25*MAC_pylon;
xc_gearpod =0.25*GP_length;


for n=1:length(machsweep)

FF_wing(n)  =(1+((0.6*(xc_wing).^-1).*(tc_wing))+(100*(tc_wing.^4))).*(1.34*(machsweep(n).^0.18)*(cos(deg2rad(Sw_wing)).^0.28));
FF_Htail(n) =(1+((0.6*(xc_Htail).^-1).*(TC_HT))+(100*(TC_HT.^4))).*(1.34*(machsweep(n).^0.18)*(cos(deg2rad(Sw_ht)).^0.28));
FF_Vtail(n) =(1+((0.6*(xc_Vtail).^-1).*(TC_VT))+(100*(TC_VT.^4))).*(1.34*(machsweep(n).^0.18)*(cos(deg2rad(Sw_vt)).^0.28));
FF_pylon(n) =(1+((0.6*(xc_pylon).^-1).*(TC_pylon))+(100*(TC_pylon.^4))).*(1.34*(machsweep(n).^0.18)*(cos(deg2rad(Sw_pylon)).^0.28));
FF_gearpod(n)=(1+((0.6*(xc_gearpod).^-1).*(TC_GP))+(100*(TC_GP.^4))).*(1.34*(machsweep(n).^0.18)*(cos(deg2rad(Sw_GP)).^0.28));
Sr_fuselage=fuse_length/Fuse_width;
Sr_nacelle=nacelle_length/nacelle_width;

FF_fuse(n)=(1+(60*Sr_fuselage.^-3)+(Sr_fuselage*400.^-1));
FF_nacelle(n)=1+(0.35*Sr_nacelle.^-1);
end



FF=[FF_wing;FF_fuse;FF_Htail; FF_Vtail ;FF_pylon;FF_nacelle ];



%interference calculation
Qc_fuse    =1.0;
Qc_wing    =1.0;
Qc_nacelle =1.3;
Qc_Htail   =1.06;
Qc_pylon   =1;
Qc_GP      =1;
Qc_Vtail   =1.06;

Qc=[Qc_fuse;Qc_wing;Qc_nacelle;Qc_Htail;Qc_Htail;Qc_GP];



CD0_wing=Cf_wing.*FF_wing.*Qc_wing.*Swet_wing.*Sref.^-1;
CD0_fuselage=Cf_fuselage.*FF_fuse.*Qc_fuse.*Swet_fuse.*Sref.^-1;
CD0_htail=Cf_HT.*FF_Htail.*Qc_Htail.*Swet_ht.*Sref.^-1;
CD0_vtail=Cf_VT.*FF_Vtail.*Qc_Vtail.*Swet_vt.*Sref.^-1;
CD0_nacelle=Cf_nacelle.*FF_nacelle.*Qc_nacelle.*Swet_nacelle.*Sref.^-1;
CD0_pylon=Cf_pylon.*FF_pylon.*Qc_pylon.*Swet_pylon.*Sref.^-1;
CD0_GP=Cf_gearpod.*FF_gearpod.*Qc_GP.*Swet_GP.*Sref.^-1;

CD0=(CD0_wing+CD0_fuselage+CD0_htail+CD0_vtail+CD0_nacelle+CD0_pylon+CD0_GP);

D_parasitic=CD0.*(0.5*rho.*V.^2*Sref);

%% parasitic pc estimatiom

MTOW=36740
W=MTOW*9.81
CL=W*(0.5*rho*V.^2*Sref).^-1

AR_wing=9.8;
Sw_wing_LE=rad2deg(tan(deg2rad(Sw_wing)) + (4./AR_wing)*((0.25).*((1-lambda_wing).*(1+lambda_wing).^-1))) ;                                           
e=1.78*(1-(0.045*AR_wing.^0.68))-0.64   ; 
CDi=(CL.^2)*(pi*e*AR_wing).^-1;
Di=CDi.*(0.5*rho*V.^2*Sref);






CD0_1=[CD0_wing,CD0_fuselage,CD0_htail,CD0_vtail,CD0_nacelle,CD0_pylon,CD0_GP];
CDexcer=0.10*sum(CD0_1);
CD0_2=[CD0_wing,CD0_fuselage,CD0_htail,CD0_vtail,CD0_nacelle,CD0_pylon,CD0_GP,CDexcer];
CD0_total=(CD0_wing+CD0_fuselage+CD0_htail+CD0_vtail+CD0_nacelle+CD0_pylon+CD0_GP+CDexcer);

CD=CD0+CDi+CDexcer;
D=0.5*rho.*V.^2*Sref.*CD;
Dtotal=D_parasitic+Di;



To=61600*2; % thrust of F 70 engine (2 nos) in newtons
T=To*(rho*(rho0.^-1))*(1-(0.25*machsweep))  


figure
%% plot for task 4
plot(machsweep,D_parasitic,'.-g',machsweep,Di,'.-b',machsweep,Dtotal,'-r')
xlabel('velocity expressed in mach ')
ylabel('drag: total,induced and parasitic(Newtons)')
ylim([0 70000])
legend('Parasitic','Induced','Total')
title('Parasitic Induced and Total drag for Fokker 70 aircraft')


%figure
%%VSp validation bar chart (requires single velocity value)
%CD0_pc=((CD0_2.*CD0_total.^-1).*100)';
%CD0_vsp=[26.21,25.30,7.61,4.23,20.51,2.66,4.39,9.09]'
%Z=[CD0_vsp,CD0_pc]
%bar(Z) %%,'Wing','Fuselage','Htail','Vtail','Nacelle','Pylon','Gearpod','excercense')

%xlabel('component number from breakup table ')
%ylabel('OPENVSP and MATLAB comparison')
%legend('OpenVSP','MATLAB')
%% plot for task 5
figure
plot(machsweep,D_parasitic,'.-g',machsweep,Di,'.-b',machsweep,Dtotal,'-r',machsweep,T,'-c')
legend('Parasitic','Induced','Total','thrust line')
xlabel('velocity expressed in mach ')
ylabel('drag: total,induced and parasitic and thrust line(Newtons)')
ylim([0 70000])
title('Parasitic Induced and Total drag  with Thrust line ')
