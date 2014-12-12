load ThicknessVariationSmoothCornea_v1.mat


figure
plot([-R0*R R0*R],[E*tau0*tau(S_Stack(:,3)).*S_Stack(:,1) E*tau0*tau(S_Stack(:,3)).*S_Stack(:,1)],'k','LineWidth',2);
xlabel('Radial Coordinate (cm)','Fontsize',24);
ylabel('Radial Tension (dynes/cm)','Fontsize',24);
set(gca,'Fontsize',24);
grid on;
%xlim([-1.2 1.2])

figure
plot([-R0*R R0*R],[AngularTension AngularTension],'k','LineWidth',2);
xlabel('Radial Coordinate (cm)','Fontsize',24);
ylabel('Angular Tension (dynes/cm)','Fontsize',24);
set(gca,'Fontsize',24);
grid on;
%xlim([-1.2 1.2])

figure
plot([-R0*R R0*R],[p p],'k','LineWidth',2);
xlabel('Radial Coordinate (cm)','Fontsize',24);
ylabel('Suction Pressure (dynes/cm^2)','Fontsize',24);
set(gca,'Fontsize',24);
grid on;
%xlim([-1.2 1.2])

figure
plot([-R0*R R0*R],[(R-S_Stack(:,3))./S_Stack(:,3) (R-S_Stack(:,3))./S_Stack(:,3)],'k','LineWidth',2);
xlabel('Radial Coordinate (cm)','Fontsize',24);
ylabel('Hoop Strain','Fontsize',24);
set(gca,'Fontsize',24);
grid on;
%xlim([-1.2 1.2])

figure
plot([-R0*R R0*R],[-E*tau0*tau(S_Stack(:,3)).*S_Stack(:,2) E*tau0*tau(S_Stack(:,3)).*S_Stack(:,2)],'k','LineWidth',2);
xlabel('Radial Coordinate (cm)','Fontsize',24);
ylabel('z (dynes/cm)','Fontsize',24);
set(gca,'Fontsize',24);
grid on;
%xlim([-1.2 1.2])

figure
plot([-R0*R(1:end-1) R0*R(1:end-1)],[RadialStrain(1:end-1) RadialStrain(1:end-1)],'k','LineWidth',2);
xlabel('Radial Coordinate (cm)','Fontsize',24);
ylabel('Radial Strain','Fontsize',24);
set(gca,'Fontsize',24);
grid on;
%xlim([-1.2 1.2])

figure
plot([-R0*R R0*R],[-R0*S_Stack(:,3) R0*S_Stack(:,3)],'--k',[-R0*R R0*R],[-R0*R R0*R],'r','LineWidth',2);
xlabel('Radial Coordinate (cm)','Fontsize',24);
ylabel('Radial Coordinate (cm)','Fontsize',24);
set(gca,'Fontsize',24);
grid on;
%xlim([-1.2 1.2])

Rs = R0*R;

save FlowVariablePaper1.mat Rs p;