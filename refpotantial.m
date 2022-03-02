function [U_p, U_n]= refpotantial (theta_p, theta_n,KokamOCVNMC, KokamNMC, KokamOCVC, KokamC, OCVcell)


%% Data fit
%  U_p = ( -4.656 + 88.669*theta_p.^2 - 401.119*theta_p.^4 + 342.909*theta_p.^6 - ...
%            462.471*theta_p.^8 + 433.434*theta_p.^10)./...
%         ( -1 + 18.933*theta_p.^2 - 79.532*theta_p.^4 + 37.311*theta_p.^6 - ...
%            73.083*theta_p.^8 + 95.96*theta_p.^10);   
% 
%  U_n= 0.7222 + 0.1387*theta_n + 0.0290*theta_n.^(1/2) - 0.0172./theta_n + ...
%         0.0019./(theta_n.^(1.5)) + 0.2808*exp(0.90-15*theta_n) - ...
%         0.7984*exp(0.4465*theta_n-0.4108);
     
%% half-OCP of the electrodes
ocv_p= KokamOCVNMC(:,2);
fp= KokamOCVNMC(:,1);
U_p= interp1(fp,ocv_p,theta_p);

ocv_n= KokamOCVC(:,2);
fn= KokamOCVC(:,1);
U_n= interp1(fn,ocv_n,theta_n);

%%  soc vs OCV
% ocv_p= KokamNMC(:,2);
% fp= KokamNMC(:,1);
% U_p= interp1(fp,ocv_p,theta_p,'linear','extrap');
% 
% ocv_n= KokamC(:,2);
% fn= KokamC(:,1);
% U_n= interp1(fn,ocv_n,theta_n,'linear','extrap');
% 

    
       
end
