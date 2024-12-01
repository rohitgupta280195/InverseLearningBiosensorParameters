function y = cost(x)

%% decision variables (free parameters)
D = x(1);                   % Diffusivity (cm2/s)
k0 = x(2);                  % Rate constant (cm/s)
alfaA = x(3);               % Cathodic charge transfer coefficient
Ae = x(4);                 % reference electrode potential (volt)

%% load experimental data file
load('paperSimulation.mat','Material','Detection')
Ipox_exp = Detection(2,1);
Ipred_exp = Detection(2,2);
Epp_exp = Detection(2,3);

%% governing equation paramters (fixed)
%D = 8.45E-06;                  % Diffusivity (cm2/s)
Cb = 1E-05;                    % Bulk concentration (mol/cm3)

%% Parameter of Butler-Volmer equation (fixed)
%k0 = 7E-4;                  % Rate constant (cm/s)
%alfaA = 0.4;                % Cathodic charge transfer coefficient
alfaC = 1-alfaA;             % Anodic charge transfer coefficient
Ef0 = Detection(2,4);         % reference electrode potential (volt)
F = 96485.332;              % Faraday's constant (C/mol)
R = 8.314;                  % Universal gas constant (J/(mol-K))
T = 298;                    % Temperature (K)

%% Experimental parameters (fixed)
EL = -0.6;                    % Voltametry starting voltage
EM = 0.25;
ER = 0.1;                     % Voltametry ending voltage
v = 0.1;                    % Scan rate (V/s)
tmax = 2*(abs(EL-ER)/v);    % Total duration of experiment (s)
%Ae = 0.1256;              % Electrode Area (cm2)

%% Simulation parameters (fixed)
dE = 1E-3;                  % Potential steps (V)
dt = dE/v;                  % Duration of each time step (s)
t = dt:dt:tmax;             % Time vector
N = 500;                    % Number of grid points
L = 6*sqrt(D*tmax);         % Span of simulation domain (cm)
dx = L/(N-1);               % Space steps (cm)
x = 0:dx:L;                 % Space vector

%% Initialization
Ci = Cb*ones(N,1);         % Initial concentraion field
CO = Ci;                   % Initialization of concentration vector reactant
CR = zeros(N,1);           % Initialization of concentration vector product
E(1) = ER;                 % Initialization of potential vector
I = zeros(1,length(t));

%% Voltage vector evaluation for transient boundary condition
for m=1:length(t)
     if m < (0.5*length(t))
        E(m+1) = E(m) - dE;
     else 
        E(m) = E(m-1) + dE;
    end
    % if m < (0.25*length(t))
    %     E(m+1) = E(m) - dE;
    % elseif m < (0.75*length(t))
    %     E(m) = E(m-1) + dE;
    % else
    %     E(m) = E(m-1) - dE;
    % end
end
% plot(t,E)
% return

%% TDMA coefficients
aP = dx/dt;               % Mass diffusion coeff for time neighbour
aE = D/dx;                 % Mass diffusion for east face
aW = D/dx;                 % Mass diffusion for west face

%% Main time loop
for m=1:length(t)

    %solver definition Mat x CO = Source

    % Lefthand matrix generation (Mat)
    % for specie A
    MatA(1,1) = aE + aP + k0*exp(-(((alfaA*F)/(R*T))*(E(m)-Ef0)));
    MatA(1,2) = -aE;
    for i = 2:N-1
        MatA(i,i) = aE + aW + aP;    % coefficient of t+dt term at P node of left side
        MatA(i,i-1) = -aW;           % coefficient of t+dt term at W node of left side
        MatA(i,i+1) = -aE;           % coefficient of t+dt term at W node of left side
    end
    MatA(N,N-1) = -aW;
    MatA(N,N) = aW + aP + ((2*D)/dx);

    %for specie B
    MatB(1,1) = aE + aP + k0*exp(((alfaC*F)/(R*T))*(E(m)-Ef0));
    MatB(1,2) = -aE;
    for i = 2:N-1
        MatB(i,i) = aE + aW + aP;    % coefficient of t+dt term at P node of left side
        MatB(i,i-1) = -aW;           % coefficient of t+dt term at W node of left side
        MatB(i,i+1) = -aE;           % coefficient of t+dt term at W node of left side
    end
    MatB(N,N-1) = -aW;
    MatB(N,N) = aW + aP + ((2*D)/dx);


    % Righthand matrix generation (Source)
    % for specie A
    SourceA(1) = aP*CO(1) + CR(1)*k0*exp(((alfaC*F)/(R*T))*(E(m)-Ef0));
    for i = 2:N-1
        SourceA(i) = aP*CO(i);    % Righthand side term at at P node
    end
    SourceA(N) = aP*CO(N) + ((2*D*Cb)/dx);
    % for specie B
    SourceB(1) = aP*CR(1) + CO(1)*k0*exp(-(((alfaA*F)/(R*T))*(E(m)-Ef0)));
    for i = 2:N-1
        SourceB(i) = aP*CR(i);    % Righthand side term at at P node
    end
    SourceB(N) = aP*CR(N) ;

    %matrix inversion and current
    CO = SourceA/MatA;
    CR = SourceB/MatB;
    I(m) = (Ae*F*D*(1/dx)*(CO(1)-CO(2)))*10^6; %current in microAmp;

    % mole fraction and surface concentration
    XO = CO(1)/Cb;
    XR = CR(1)/Cb;
    XO_surface(m) = CO(1)/Cb;
    XR_surface(m) = CR(1)/Cb;

%     if ((XO_surface(m)+XR_surface(m))>1)
%         disp('Mass Conservation Fails')
%     end

%     figure(1)
%     set(gcf,'WindowState','maximized')
%     set(gcf,'color','w');
%     tiledlayout(3,2,'TileSpacing','tight','Padding','loose')
%     
%     nexttile
%     plot(t,E,'k','LineWidth',1.5)
%     hold on
%     plot(t(m),E(m),'or','MarkerFaceColor','r')
%     hold off
%     grid on
%     xlabel('{\it t} (sec)','FontWeight','bold')
%     ylabel('{\it E} (volt)','FontWeight','bold')
%     xlim([0 tmax])
%     ylim([EL ER])
%     grid on
% 
%     nexttile
%     plot(t(1:m),I(1:m),'k','LineWidth',1.5)
%     hold on
%     plot(t(m),I(m),'or','MarkerFaceColor','r')
%     hold off
%     grid on
%     xlabel('{\it t} (sec)','FontWeight','bold')
%     ylabel('{\it I} (\muA)','FontWeight','bold')
%     xlim([0 tmax])
%     ylim([-30 +30])
% 
%     nexttile
%     grid on
%     plot((x*10),CO,'r','LineWidth',1.5);
%     hold on
%     plot((x*10),CR,'b','LineWidth',1.5);
%     xlabel('{\it Z} (mm)','FontWeight','bold')
%     ylabel('{\it C} (mol/cm^{3})','FontWeight','bold')
%     legend('[Fe(CN)_6]^{3-}','[Fe(CN)_6]^{4-}')
%     xlim([0 (L*10)])
%     %ylim([0 Cb])
%     hold off
%     grid on
% 
%     nexttile
%     plot(E(1:m),I(1:m),'k','LineWidth',1.5)
%     hold on
%     plot(E(m),I(m),'or','MarkerFaceColor','r')
%     hold off
%     xlim([EL ER])
%     ylim([-30 +30])
%     xlabel('{\it E} (volt)','FontWeight','bold')
%     ylabel('{\it I} (\muA)','FontWeight','bold')
%     grid on
% 
%     nexttile
%     Cat = categorical({'[Fe(CN)_6]^{3-}','[Fe(CN)_6]^{4-}'});
%     Y = [XO XR];
%     bar(Cat,Y);
%     xlabel('Species','FontWeight','bold')
%     ylabel('Mole Frac','FontWeight','bold')
%     ylim([0 1])
%     grid on
% 
%     nexttile
%     plot(t(1:m),XO_surface(1:m),'r','LineWidth',1.5)
%     hold on
%     plot(t(1:m),XR_surface(1:m),'b','LineWidth',1.5)
%     xlabel('{\it t} (sec)','FontWeight','bold')
%     ylabel('Mole Frac','FontWeight','bold')
%     %legend('[Fe(CN)_6]^{3-}','[Fe(CN)_6]^{4-}')
%     grid on
%     %ylim([0 1])
%     xlim([0 tmax])
%     drawnow
%     Frame(m) = getframe(gcf);
end

% 
% % create the video writer with 1 fps
%   writerObj = VideoWriter('myVideo.avi');
%   writerObj.FrameRate = 10;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(Frame)
%     % convert the image to a frame
%     frame = Frame(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);



    % figure(1)
    % plot(E,I,'k','LineWidth',1.5)
    % %xlim([EL ER])
    % %ylim([-100 +100])
    % xlabel('{\it E} (volt)','FontWeight','bold')
    % ylabel('{\it I} (\muA)','FontWeight','bold')
    % grid on

[Ipox_sim,indMaxI] = max(I);
[Ipred_sim,indMinI] = min(I);
Epp_sim = (E(indMaxI) - E(indMinI));

y(1) = abs(((Ipox_sim - Ipox_exp)/Ipox_exp)*100);
y(2) = abs(((Ipred_sim - Ipred_exp)/Ipred_exp)*100);
y(3) = abs(((Epp_sim - Epp_exp)/Epp_exp)*100);


end