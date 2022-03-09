clc
clear all; 

%% setting
num = 1*1e4;          % simulation number             
P_dBM = 15:25;              %transmit power dBm
P_t= 10.^((P_dBM-30)./10);

BW=10*10^6; %10 MHz
Nf=10;%dB
sigma2_dbm= -170+10*log10(BW)+Nf; %Thermal noise in dBm
sigma_square=10^((sigma2_dbm-30)/10);

fc = 1e9; % 1 GHz carrier
c = 3*10^8; % speed of light
wavelength=3*10^8./fc;%wavelength
C_L=(wavelength/(4*pi))^2; %intercept of NLOS


rho = P_t./sigma_square; %%transmit power
rho_dB = 10.*log10(rho);

alpha = 2.4;                 % large-scale parameter

N=30; % num of RIS elements

k = 2; %% Rician parameter
s = sqrt(k./(1+k)); % matlab paramiter one
sigma = sqrt( 1./ (2.*(1+k))  ); % matlab paramiter Two

a_rfr = 0.6;
a_rfl = 0.4; %power allocation from BS

beta_rfr = 0.3; %power allocation by RIS
beta_rfl = 0.7;

R = 20; % RIS surving range
r_1 = 100; % distance between BS to RIS

R_user=0.1;  % threshold of SNR
R_SIC=0.1;
gamma_th =2.^(R_user)-1;
gamma_th_SIC=2.^(R_SIC)-1;
% 
A= 30;
B_rfr = 22.46;% 6.739;
B_rfl = 22.46;% 15.58 ;


% A= 20;
% B_rfr = 15.06;
% B_rfl = B_rfr;

Gamma = max([gamma_th_SIC*sigma_square/(a_rfr-gamma_th_SIC*a_rfl), gamma_th*sigma_square/a_rfl]);
Gamma3 = gamma_th*sigma_square/(a_rfr-gamma_th*a_rfl);
% %% simulation 
% P_out_sim_rfr = zeros(1,length(rho));
% P_out_sim_rfl = zeros(1,length(rho));
% for i = 1:length(rho)
%     tic;
%     r_rfr = sqrt(R^2.*rand(1,num));          % user position :reflection
%     r_rfl = sqrt(R^2.*rand(1,num));          % user position :refraction
%     Poutsum_rfr = 0;
%     Poutsum_rfl = 0;
%     for j =1:1:num
%     g_rfr = beta_rfr.*random('Gamma',A,B_rfr,[1,1]);
%     g_rfl = beta_rfl.*random('Gamma',A,B_rfl,[1,1]);
%     
%     SNR_SIC(j) = a_rfr*P_t(i)*C_L*r_1^(-alpha)*C_L*r_rfl(j)^(-alpha)*g_rfl / ...
%                  (a_rfl*P_t(i)*C_L*r_1^(-alpha)*C_L*r_rfl(j)^(-alpha)*g_rfl +sigma_square);
%     SNR_rfl(j) = a_rfl*P_t(i)*C_L*r_1^(-alpha)*C_L*r_rfl(j)^(-alpha)*g_rfl /sigma_square ;
%     SNR_rfr(j) = a_rfr*P_t(i)*C_L*r_1^(-alpha)*C_L*r_rfr(j)^(-alpha)*g_rfr / ...
%                  (a_rfl*P_t(i)*C_L*r_1^(-alpha)*C_L*r_rfr(j)^(-alpha)*g_rfr +sigma_square);
%              
%          if   SNR_SIC(j) < gamma_th_SIC   || SNR_rfl(j) < gamma_th 
%              Poutsum_rfl= Poutsum_rfl+1;
%          end
%          if SNR_rfr(j) < gamma_th
%              Poutsum_rfr = Poutsum_rfr +1;
%          end
%      
%     end
%      P_out_sim_rfr(i) = Poutsum_rfr/num;
%      P_out_sim_rfl(i) = Poutsum_rfl/num;
%     
%     toc;
% end
% 
% figure
% semilogy(rho_dB,P_out_sim_rfr,'bo-','LineWidth',2);
% hold on;
% 
% semilogy(rho_dB,P_out_sim_rfl,'r*-','LineWidth',2);
% hold on;

%% analysis with curve fitting
for i = 1:length(rho)
    if a_rfr >gamma_th_SIC*a_rfl && a_rfr >gamma_th*a_rfl
        index_rfr = Gamma3/(P_t(i)*C_L^2)*r_1^alpha/B_rfr/beta_rfr ;
        index_rfl = Gamma/(P_t(i)*C_L^2)*r_1^alpha/B_rfl/beta_rfl ;
        P_out_rfl(i) = 2/R^2*integral(@(x)x.*gammainc(index_rfl.*x.^alpha,A) ,0,R);
        P_out_rfr(i) = 2/R^2*integral(@(x)x.*gammainc(index_rfr.*x.^alpha,A) ,0,R);
    end
end

semilogy(rho_dB,P_out_rfl,'r-')
hold on;
semilogy(rho_dB,P_out_rfr,'b-')
hold on;
