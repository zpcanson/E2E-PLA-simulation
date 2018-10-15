clear all
clc
close all

%% (Signal Source) Generate a random binary data stream
N = 100000;   % number of bits


addpath('E:\program test\channels\functions');
%% the number of path
L1=2;
L2=3;
%L=L1*L2;
L=6;
deta_h_Value=0.2:0.2:3;
deta_d_Value=0.2:0.2:3;

% deta_h_Value=2:-0.2:0;
% deta_d_Value=2:-0.2:0;

Z_Value = 1;
kappa_h_Value= 0 ;
kappa_tau_Value = 0;
%% scenarios  
    fade_scn=1;
    SNR_dB=[10];% SNR scan  
    SNR=10.^(SNR_dB/10);
   
%% distance between two channel uses
    ch_dis=1; % 
    ch_scn=1;
% normalized Dopplers frequencies


vfar = [.001,  .02,  .05 ];
vfrb = [.001,  .03,  .05 ];

vfd_ar=[.001,  .02,  .05 ];
vfd_rb=[.001,  .03,  .05 ];

vfer = [.001,  .001,  .001 ];
vfd_er=[.001,  .001,  .001 ];
% scenario number
fer=vfer(fade_scn);
far=vfar(fade_scn);
frb=vfrb(fade_scn);

fd_ar=vfd_ar(fade_scn);
fd_rb=vfd_rb(fade_scn);
fd_er=vfd_er(fade_scn);

% auto-correlation values

alfa1=besselj(0,2*pi*far*ch_dis);
alfa2=besselj(0,2*pi*frb*ch_dis);
alfa=alfa1*alfa2;

rho1=besselj(0,2*pi*fd_ar*ch_dis);
rho2=besselj(0,2*pi*fd_rb*ch_dis);

Pmd__SNR8 = zeros(length(deta_h_Value),length(deta_d_Value),length(SNR_dB));
Pfa_SNR8  = zeros(length(deta_h_Value),length(deta_d_Value),length(SNR_dB));


for i=1:length(SNR_dB)
    
    P_theoretical = zeros(length(deta_h_Value),length(deta_d_Value));
    Pmd_theory_h_d = zeros(length(deta_h_Value),length(deta_d_Value));
    Pfa_theoretical = zeros(length(deta_h_Value),length(deta_d_Value));

    Pfa_theoretical = zeros(length(deta_h_Value),length(deta_d_Value));

 for h_indx=1:length(deta_h_Value)
    Z=1;
    kappa_h = 0;
    kappa_tau = 0;
    deta_h = deta_h_Value(h_indx);
   
    
 for d_indx=1:length(deta_d_Value)
     
   deta_d = deta_d_Value(d_indx);
   O_d=0;
   O_h=0;
   S_h=0;
   S_d=0;
   S=0;
   Pd_theoret=0;
   numerr=0;
   
   sig_ar2=1;    %define SNR=sigma^2/sig_w2 
   sig_rb2=1;    %define SNR=sigma^2/sig_w2 
   sig_er2=sig_ar2*10.^(kappa_h/10);    %define SNR=sigma^2/sig_w2 
    %define SNR=sigma^2/sig_w2 
   sig_w2 =1/SNR(i);       % noise variance
   
    hrd=sqrt(sig_rb2).*flat_cos(N,frb,ch_dis);
    h2=abs(hrd).^2;
    h2_mean=mean(h2);
    pow_h2= mean(hrd.*conj(hrd));
    ear=sqrt(sig_ar2)*cxn(N,sig_w2 );
    ear2=abs(ear).^2;
    ear2_mean=mean(ear2);
    pow_h2= mean(hrd.*conj(hrd));
 
    Ph_H11=exp(-deta_h./((alfa2^2*sig_er2+sig_ar2).*h2+(1-alfa2^2)*sig_er2.*sig_rb2+2*sig_w2 ));
    Ph_H1 =mean( Ph_H11);
    
%   Ph_H2=exp(-deta_h./((alfa2^2*sig_er2+sig_ar2).*sig_rb2+(1-alfa2^2)*sig_er2*sig_rb2+2*sig_w2 ));
%   Ph_H33=exp(-deta_h./((alfa2^2*sig_er2+sig_ar2).*h2+(1-alfa2^2)*sig_er2.*ear2+2*sig_w2 ));
%   Ph_H3=mean( Ph_H33); 

   sig_tau_ar2=1;   %define SNR=sigma^2/sig_w^2
   sig_tau_rb2=0.8;   %define SNR=sigma^2/sig_w^2
   
  % sig_tau_er2=sig_tau_ar2*10.^(kappa_tau/10);   %define SNR=sigma^2/sig_w^2
   sig_tau_er2=sig_tau_ar2*10.^(kappa_tau/10);  
   s1 =1/(sig_tau_er2+ sig_tau_ar2+sig_w2 );
   s2 =1/sqrt(4*(sig_tau_rb2)^2*(1-rho2^2) + 4*(sig_tau_rb2)*sig_w2 +sig_w2 ^2);
   F =1+s1 *s2 /(s1 .^2-s2 .^2)*(s2 /s1 *exp(-s1 *deta_d)-...
         s1 /s2 *exp(-s2 *deta_d));
   
   P_tau_H1 =1-F ;

    Sigma_EA = sig_tau_er2 +  sig_tau_ar2 + sig_w2 ;
    Sigma_RB = sqrt(4*sig_tau_rb2^2*(1-rho2^2) + 4*sig_tau_rb2* sig_w2 + sig_w2 ^2 );
 %   F_tau_H1 =  Sigma_RB/Sigma_EA/(1-(Sigma_RB/Sigma_EA)^2)*(Sigma_EA/Sigma_RB*(1-exp(-Sigma_RB*deta_d)) -Sigma_RB/Sigma_EA*(1-exp(-Sigma_EA*deta_d)) );
%    P_tau_H1 = 1/(1-(Sigma_EA/Sigma_RB)^2)*((Sigma_EA/Sigma_RB)^2*exp(- deta_h/Sigma_EA)-exp(- deta_h/Sigma_RB));
%
    F_delt_EB  = 1/(1-(Sigma_RB/Sigma_EA)^2)*((1-exp(- Sigma_RB*deta_h)- (Sigma_RB/Sigma_EA)^2*(1-exp(- Sigma_EA*deta_h))));
    P_tau_H21  = 1- F_delt_EB ;
    
    S1 =sig_tau_er2 +  sig_tau_ar2 + sig_w2 ;
    S2 =sqrt(4*(sig_tau_rb2)^2*(1-rho2^2) + 4*(sig_tau_rb2)*sig_w2 +sig_w2 ^2);
    F_S1S2 = S2 /S1 *1/(1-(S2 /S1 )^2)*(S1 /S2 *(1-exp(- S2 *deta_d))-S2 /S1 *(1-exp(- S1 *deta_d)));
    P_tau_H11  = 1- F_S1S2 ;
    
  if Z < L-1 
   for z=0:Z
      for v=0:z
      P_theoretical(h_indx,d_indx)=P_theoretical(h_indx,d_indx)+nchoosek(L,v)*Ph_H1 ^(v)*(1-Ph_H1 )^(L-v)*nchoosek(L-1,z-v)*(P_tau_H1 )^(z-v)*(1-P_tau_H1 )^(L-1-z+v);
      Pmd_theory_h_d(h_indx,d_indx)= P_theoretical(h_indx,d_indx);
      end
   end     
  end
     
  if Z>=L-1  
    for z=0:L-1
      for v=0:z
      P_theoretical(h_indx,d_indx) = P_theoretical(h_indx,d_indx)+nchoosek(L,v)*Ph_H1 ^(v)*(1-Ph_H1 )^(L-v)*nchoosek(L-1,z-v)*P_tau_H1 ^(z-v)*(1-P_tau_H1 )^(L-1-z+v);
      end
    end
    
    for z=L:Z
      for v=z-L:L-1
      P_theoretical(h_indx,d_indx) = P_theoretical(h_indx,d_indx) + nchoosek(L,z-v)*Ph_H1 ^(z-v)*(1-Ph_H1 )^(L-z+v)*nchoosek(L-1,v)*(P_tau_H1 )^(v)*(1-P_tau_H1 )^(L-1-v);
   %   Pmd_theory_h_d(ka_indx,i) = P_theoretical(ka_indx,i);
      end
    end   
  end

      
      Pmd__SNR8(h_indx,d_indx,i)=P_theoretical(h_indx,d_indx);

      
%  for m=1:L1*L2
%    for j=1:L1
%        har(j,:)=sqrt(sig_ar2).*flat_cos(N,far,ch_dis);
%        her(j,:)=sqrt(sig_er2).*flat_cos(N,fer,ch_dis);
% 
%      for z=1:L2
%        hrb(z,:)=sqrt(sig_rb2).*flat_cos(N,frb,ch_dis);
%        harb(m,:)=har(j,:).*hrb(z,:);
%   
%        herb(m,:)=her(j,:).*hrb(z,:);
% 
%        z_ar(m,:)=cxn(N,sig_w2 );
%        z_er(m,:)=cxn(N,sig_w2 );
%        
%        harb_hat(m,:)=harb(m,:) + z_ar(m,:); 
%        herb_hat(m,:)=herb(m,:) + z_er(m,:);
%      end
%    end
%  end
%            
%    for ii=1:L-1   
%         % AWGN noise CN(0,sig_w2 )
%         z_dR_ar(ii,:)=sqrt(sig_w2 /2)*randn(1,N);
%         z_dI_ar(ii,:)=sqrt(sig_w2 /2)*randn(1,N);
%         
%         z_dR_rb(ii,:)=sqrt(sig_w2 /2)*randn(1,N);
%         z_dI_rb(ii,:)=sqrt(sig_w2 /2)*randn(1,N);
%         
%         z_dR_er(ii,:)=sqrt(sig_w2 /2)*randn(1,N);
%         z_dI_er(ii,:)=sqrt(sig_w2 /2)*randn(1,N);
% 
%         power_zR_a(ii)=mean(z_dR_ar(ii,:).*conj(z_dR_ar(ii,:)));
%         power_zI_a(ii)=mean(z_dI_ar(ii,:).*conj(z_dI_ar(ii,:)));
%  
%         % generate time delay    
%         di_ar(ii,:)=sqrt(2*sig_tau_ar2)*flat_cos(N,fd_ar,ch_dis);
%         di_rb(ii,:)=sqrt(2*sig_tau_rb2)*flat_cos(N,fd_rb,ch_dis);
% 
%         di_er(ii,:)=sqrt(2*sig_tau_er2)*flat_cos(N,fd_er,ch_dis);
%         
%         power_di(ii)=mean(di_ar(ii,:).*conj(di_ar(ii,:)));
% 
%         d_R_ar(ii,:) = real(di_ar(ii,:));
%         d_I_ar(ii,:) = imag(di_ar(ii,:));
%         
%         d_R_rb(ii,:) = real(di_rb(ii,:));
%         d_I_rb(ii,:) = imag(di_rb(ii,:));
% 
%         d_R_er(ii,:)=real(di_er(ii,:));
%         d_I_er(ii,:)=imag(di_er(ii,:));
%         
%         power_d_R_a(ii)=mean(d_R_ar(ii,:).*conj(d_R_ar(ii,:)));
%         power_d_I_a(ii)=mean(d_I_ar(ii,:).*conj(d_I_ar(ii,:)));
% 
%         d_a(ii,:)=d_R_ar(ii,:).^2 + d_I_ar(ii,:).^2;
%         
%         % estimation delay
%         d_R_ar_hat(ii,:)=d_R_ar(ii,:) + z_dR_ar(ii,:);
%         d_I_ar_hat(ii,:)=d_I_ar(ii,:) + z_dI_ar(ii,:);
%         
%         d_R_rb_hat(ii,:)=d_R_rb(ii,:) + z_dR_rb(ii,:);
%         d_I_rb_hat(ii,:)=d_I_rb(ii,:) + z_dI_rb(ii,:);
%         
%         d_R_er_hat(ii,:)=d_R_er(ii,:) + z_dR_er(ii,:);
%         d_I_er_hat(ii,:)=d_I_er(ii,:) + z_dI_er(ii,:);
%         
%         power_d_R_a_hat(ii)=mean(d_R_ar_hat(ii,:).*conj(d_R_ar_hat(ii,:)));
%         power_d_I_a_hat(ii)=mean(d_I_ar_hat(ii,:).*conj(d_I_ar_hat(ii,:)));
% 
%         d_ar_hat(ii,:)=d_R_ar_hat(ii,:).^2 + d_I_ar_hat(ii,:).^2;
%         d_rb_hat(ii,:)=d_R_rb_hat(ii,:).^2 + d_I_rb_hat(ii,:).^2;
%         
%         d_er_hat(ii,:) = d_R_er_hat(ii,:).^2 + d_I_er_hat(ii,:).^2;
%         
%     end    
%       
%  for n=2:N     
%      for j=1:L
%         Q_h(j,n-1)=abs(herb_hat(j,n)-harb_hat(j,n-1)).^2;
%          if  Q_h(j,n-1) > deta_h;
%              O_h(j,n-1)=1;
%          else O_h(j,n-1)=0;
%          end  
%      end
%      S_h(n-1)=sum(O_h(:,n-1));
% 
%      
%     for ii=1:L-1
%       Q_d(ii,n-1)=abs(d_er_hat(ii,n) + d_rb_hat(ii,n) - d_ar_hat(ii,n-1) - d_rb_hat(ii,n-1));       
%          if  Q_d(ii,n-1) > deta_d;
%              O_d(ii,n-1)=1;
%          else O_d(ii,n-1)=0;
%          end
%     end
%          
%     S_d(n-1)=sum(O_d(:,n-1));
%     S(n-1)=S_h(n-1)+S_d(n-1);
% 
%     if S(n-1)> Z
%          numerr(i,n-1)=1;
%     else numerr(i,n-1)=0;
%     end
%      
%  end 
%     
%  % Pfa_theoretical =Pfa_theoretical1;
%   num_pd(ka_indx,i)=sum(numerr(i,:));
%   Pmd_simu_h_d(ka_indx,i)=1-num_pd(ka_indx,i)/(N-1);
end

%% scenarios  
% auto-correlation values

alfa1=besselj(0,2*pi*far*ch_dis);
alfa2=besselj(0,2*pi*frb*ch_dis);
alfa=alfa1*alfa2;
a=alfa;
rho1=besselj(0,2*pi*fd_ar*ch_dis);
rho2=besselj(0,2*pi*fd_rb*ch_dis);


for d_indx=1:length(deta_d_Value)
     
   deta_d = deta_d_Value(d_indx);


   O_h=0;
   O_d=0;
   numerr=0;

    % hrd=sqrt(sig_rb )*cxn(N,sig_w2 );
    hrd=sqrt(sig_rb2).*flat_cos(N,frb,ch_dis);
    h2=abs(hrd).^2;
    h2_mean=mean(h2);
    pow_h2= mean(hrd.*conj(hrd));
    ear=sqrt(sig_ar2)*cxn(N,sig_w2 );
    ear2=abs(ear).^2;
    ear2_mean=mean(ear2);
    pow_h2= mean(hrd.*conj(hrd));
    Pb_H01=exp(-deta_h./(2*(1-alfa)*sig_ar2.*h2+2*sig_w2 ));
    Pb_H02=exp(-deta_h./(2*(1-alfa)*sig_ar2.*sig_rb2+2*sig_w2 ));
    Pb_H03=exp(-deta_h./((1-a)^2*sig_ar2.*sig_rb2+(1-a^2).*h2.*ear2+2*sig_w2 ));
    Pb_H04=exp(-deta_h./((1-a)^2*sig_ar2.*sig_rb2+(1-a^2).*sig_rb2.*ear2+2*sig_w2 ));
 %   Pb_H05=exp(-deta_h./((1-a)^2*sig_ar2.*h2+(1-a^2).*h2.*ear2+2*sig_w2 ));
 %   Pb_H06=exp(-deta_h./((1-a)^2*sig_ar2.*h2+(1-a^2).*h2.*sig_ar2+2*sig_w2 ));

    Ph_H0 = mean(Pb_H02);
    
 %  Ph_H0 =exp(-deta_h/(2*(1-alfa1)*sig_a .^2+2*sig_w^2));
    Pd_H0 =exp(-deta_d/sqrt(4*(sig_tau_ar2)^2*(1-rho1^2) + 4*(sig_tau_ar2)*sig_w2 +sig_w2 ^4));

   s1 =1/sqrt(4*(sig_tau_ar2)^2*(1-rho1^2) + 4*(sig_tau_ar2)*sig_w2 +sig_w2 ^2);
   s2 =1/sqrt(4*(sig_tau_rb2)^2*(1-rho2^2) + 4*(sig_tau_rb2)*sig_w2 +sig_w2 ^2);
   F =1+s1 *s2 /(s1 .^2-s2 .^2)*(s2 /s1 *exp(-s1 *deta_d)-...
         s1 /s2 *exp(-s2 *deta_d));
   
   P_tau_H0 =1-F ;    
    
    Sigma_AR = sqrt(4*sig_tau_ar2^2*(1-rho1^2) + 4*sig_tau_ar2* sig_w2 + sig_w2 ^2 );
    Sigma_RB = sqrt(4*sig_tau_rb2^2*(1-rho2^2) + 4*sig_tau_rb2* sig_w2 + sig_w2 ^2 );
 %   F_tau_H1 =  Sigma_RB/Sigma_EA/(1-(Sigma_RB/Sigma_EA)^2)*(Sigma_EA/Sigma_RB*(1-exp(-Sigma_RB*deta_d)) -Sigma_RB/Sigma_EA*(1-exp(-Sigma_EA*deta_d)) );
%    P_tau_H1 = 1/(1-(Sigma_EA/Sigma_RB)^2)*((Sigma_EA/Sigma_RB)^2*exp(- deta_h/Sigma_EA)-exp(- deta_h/Sigma_RB));
%
    F_delt_AB  = 1/(1-(Sigma_RB/Sigma_AR)^2)*((1-exp(- Sigma_RB*deta_d)- (Sigma_RB/Sigma_AR)^2*(1-exp(- Sigma_AR*deta_d))));
    P_tau_H000  = 1- F_delt_AB ;
   
   
    
    Sar =sqrt(4*(sig_tau_ar2)^2*(1-rho1^2) + 4*(sig_tau_ar2)*sig_w2 +sig_w2 ^2);
    Srb =sqrt(4*(sig_tau_rb2)^2*(1-rho2^2) + 4*(sig_tau_rb2)*sig_w2 +sig_w2 ^2);
    F_Sar_Srb = Srb /Sar *1/(1-(Srb /Sar )^2)*(Sar /Srb *(1-exp(- Srb *deta_d))-Srb /Sar *(1-exp(- Sar *deta_d)));
    P_tau_H00  = 1- F_Sar_Srb ;
    
    
   
   
   Q_h=0;
   O_h=0;
   S_h=0;
   
  if Z<=L-1 
   for z=Z+1:L-1
      for v=0:z
      Pfa_theoretical(h_indx,d_indx)=Pfa_theoretical(h_indx,d_indx)+nchoosek(L,v)*Ph_H0 ^(v)*(1-Ph_H0 )^(L-v)*nchoosek(L-1,z-v)*(P_tau_H0 )^(z-v)*(1-P_tau_H0 )^(L-1-z+v);
      end
   end
   
     for z=L:2*L-1
      for v=z-L:L-1
      Pfa_theoretical(h_indx,d_indx)=Pfa_theoretical(h_indx,d_indx)+nchoosek(L,z-v)*Ph_H0 ^(z-v)*(1-Ph_H0 )^(L-z+v)*nchoosek(L-1,v)*(P_tau_H0 )^(v)*(1-P_tau_H0 )^(L-1-v);
      end
     end
  end
     
  if Z>L-1  
    for z=Z+1:2*L-1
      for v=z-L:L-1
      Pfa_theoretical(h_indx,d_indx)=Pfa_theoretical(h_indx,d_indx)+nchoosek(L,z-v)*Ph_H0 ^(z-v)*(1-Ph_H0 )^(L-z+v)*nchoosek(L-1,v)*P_tau_H0 ^(v)*(1-P_tau_H0 )^(L-1-v);
      end
    end
  end
  Pfa_SNR8(h_indx,d_indx,i)= Pfa_theoretical(h_indx,d_indx);
  b=1;
%    for m=1:L1*L2
%    for j=1:L1
%        har(j,:)=sqrt(sig_ar2).*flat_cos(N,far,ch_dis);
%      for z=1:L2
%        hrb(z,:)=sqrt(sig_rb2).*flat_cos(N,frb,ch_dis);
%        harb(m,:)=har(j,:).*hrb(z,:);
%        z_ar(m,:)=cxn(N,sig_w2 );
%        harb_hat(m,:)=harb(m,:) + z_ar(m,:);
%      end
%    end
%  end
%  
%     
%       for ii=1:L-1
%         % AWGN noise CN(0,sig_w2 )
%         z_dR_ar(ii,:) = sqrt(sig_w2 /2)*randn(1,N);
%         z_dI_ar(ii,:) = sqrt(sig_w2 /2)*randn(1,N);
%         
%         z_dR_rb(ii,:) = sqrt(sig_w2 /2)*randn(1,N);
%         z_dI_rb(ii,:) = sqrt(sig_w2 /2)*randn(1,N);
% 
% %         power_zR_a(ii)=mean(z_dR_ar(ii,:).*conj(z_dR_ar(ii,:)));
% %         power_zI_a(ii)=mean(z_dI_a(ii,:).*conj(z_dI_a(ii,:)));
% %  
%         % generate time delay    
% 
%         di_ar(ii,:)=sqrt(2*sig_tau_ar2)*flat_cos(N,fd_ar,ch_dis);
%         di_rb(ii,:)=sqrt(2*sig_tau_rb2)*flat_cos(N,fd_rb,ch_dis);
% %         power_di(ii)=mean(di_ar(ii,:).*conj(di_ar(ii,:)));
% 
%         d_R_ar(ii,:) = real(di_ar(ii,:));
%         d_I_ar(ii,:) = imag(di_ar(ii,:));
%         
%         d_R_rb(ii,:) = real(di_rb(ii,:));
%         d_I_rb(ii,:) = imag(di_rb(ii,:));
%         
% %         power_d_R_a(ii) = mean(d_R_ar(ii,:).*conj(d_R_ar(ii,:)));
% %         power_d_I_a(ii) = mean(d_I_ar(ii,:).*conj(d_I_ar(ii,:)));
% 
%         d_ar(ii,:) = d_R_ar(ii,:).^2 + d_I_ar(ii,:).^2;
%         d_rb(ii,:) = d_R_rb(ii,:).^2 + d_I_rb(ii,:).^2;
%         
%         % estimation delay
%         d_R_ar_hat(ii,:) = d_R_ar(ii,:) + z_dR_ar(ii,:);
%         d_I_ar_hat(ii,:) = d_I_ar(ii,:) + z_dI_ar(ii,:);
%         
%         d_R_rb_hat(ii,:) = d_R_rb(ii,:) + z_dR_rb(ii,:);
%         d_I_rb_hat(ii,:) = d_I_rb(ii,:) + z_dI_rb(ii,:);
%         
%         d_ar_hat(ii,:) = d_R_ar_hat(ii,:).^2 + d_I_ar_hat(ii,:).^2;
%         d_rb_hat(ii,:) = d_R_rb_hat(ii,:).^2 + d_I_rb_hat(ii,:).^2;     
%       end
  
%  for n=2:N     
%     for j=1:L 
%         Q_h(j,n-1)=abs(harb_hat(j,n)-harb_hat(j,n-1)).^2;
%         
%          if  Q_h(j,n-1) > deta_h;
%              O_h(j,n-1)=1;
%          else O_h(j,n-1)=0;
%          end
%   
%     end
%      
%     for ii=1:L-1
%        
%         Q_d(ii,n-1)=abs(d_ar_hat(ii,n)-d_ar_hat(ii,n-1)+d_rb_hat(ii,n)-d_rb_hat(ii,n-1));
%         
%          if  Q_d(ii,n-1) > deta_d;
%              O_d(ii,n-1)=1;
%          else O_d(ii,n-1)=0;
%          end
%     end
%        
%     S_h(n-1)=sum(O_h(:,n-1));
%     S_d(n-1)=sum(O_d(:,n-1));
%     
%     S(n-1)=S_h(n-1) + S_d(n-1);
%     
%     if S(n-1)> Z
%          numerr(i,n-1)=1;
%     else numerr(i,n-1)=0;
%     end
     
  end 
    
 end
end 
 % Pfa_theoretical =Pfa_theoretical1;
%   num_fa(ka_indx,i)=sum(numerr(i,:));
%   Pfa_simu(ka_indx,i)=num_fa(ka_indx,i)/(N-1);



%plot(SNR_dB,Pfa_theoretical,'-r',SNR_dB,Pfa_simu,'ko');

% figure(1) 
% hold on; box on;
% set(0,'defaultfigurecolor','w');
% 
%   plot(SNR_dB,Pm_d(1,:),'-hk','LineWidth',2,'MarkerSize',10);
%   plot(SNR_dB,Pm_d(2,:),'pk--','LineWidth',2,'MarkerSize',10);
%   plot(SNR_dB,Pm_d(3,:),'-.sk','LineWidth',2,'MarkerSize',10);
% 
%  
%  plot(SNR_dB,Pfa_theoretical(1,:),'-+k','LineWidth',2,'MarkerSize',10);
%  plot(SNR_dB,Pfa_theoretical(2,:),'--<k','LineWidth',2,'MarkerSize',10);
%  plot(SNR_dB,Pfa_theoretical(3,:),':dk','LineWidth',2,'MarkerSize',10);
% 
%  
% %  
% %   plot(SNR_dB,Pmd_simu_h_d(1,:),'ok','LineWidth',2,'MarkerSize',10);
% %   plot(SNR_dB,Pmd_simu_h_d(2,:),'ok','LineWidth',2,'MarkerSize',10);
% %   plot(SNR_dB,Pmd_simu_h_d(3,:),'ok','LineWidth',2,'MarkerSize',10);
% %   
% %  plot(SNR_dB,Pfa_simu(1,:),'ok','LineWidth',2,'MarkerSize',10);
% %  plot(SNR_dB,Pfa_simu(2,:),'ok','LineWidth',2,'MarkerSize',10);
% %  plot(SNR_dB,Pfa_simu(3,:),'ok','LineWidth',2,'MarkerSize',10);
% % 
%  % plot(SNR_dB,Pfa_theoretical(2,:),'--<r','LineWidth',2,'MarkerSize',10);
% % plot(SNR_dB,Pfa_theoretical(3,:),':dr','LineWidth',2,'MarkerSize',10);
% 
% 
% xlabel('Average SNR per Hop [dB]','fontsize',24);
% ylabel('P_{md}/P_{fa}','fontsize',24);
 %set(gca,'XTick',SNR_dB(1):4:SNR_dB(end),'YTick',0:0.1:1,'FontSize',24,...
   % 'FontName','Times New Roman');
% legend('Z=1 P_{md}','Z=3 P_{md}','Z=5 P_{md}','Z=1 P_{fa}','Z=3 P_{fa}','Z=5 P_{fa}');
% grid on
% 
figure(1) 
[x,y]=meshgrid(deta_h_Value,deta_d_Value);
surf(x,y,Pmd__SNR8(:,:,1));
 xlabel('\delta_h','fontsize',30,'FontWeight','bold');
 ylabel('\delta_{\tau}','fontsize',30,'FontWeight','bold');
 zlabel('P_{md}','fontsize',30,'FontWeight','bold');
set(gca,'XTick',0:0.5:3,'YTick',0:0.5:3,'ZTick',0:0.05:0.3,'FontSize',24,...
   'FontName','Times New Roman','FontWeight','bold');
  
  colormap([0.5 0.5 0.5]);
  view(51,26);
  set(gca,'linewidth',1.1);
  grid on;
  set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5);
 line([0,0],[0,3],[0.3,0.3],'color','k','LineStyle','-.','linewidth',1.3)
 line([0,3],[3,3],[0.3,0.3],'color','k','LineStyle','-.','linewidth',1.3)
 line([3,3],[3,3],[0,0.3],'color','k','LineStyle','-.','linewidth',1.3)

  
  
  
  
figure(2) 
[x,y]=meshgrid(deta_h_Value,deta_d_Value);
surf(x,y,Pfa_SNR8(:,:,1));
 xlabel('\delta_h','fontsize',30,'FontWeight','bold');
 ylabel('\delta_{\tau}','fontsize',30,'FontWeight','bold');
 zlabel('P_{fa}','fontsize',30,'FontWeight','bold');
set(gca,'XTick',0:0.5:3,'YTick',0:0.5:3,'ZTick',0:0.2:1,'FontSize',24,...
   'FontName','Times New Roman','FontWeight','bold');
 view(51,26);
% set(gca,'ZTick',0:0.2:1,'FontSize',24,...
%     'FontName','Times New Roman');
 colormap([0.5 0.5 0.5]);

 grid on;
 set(gca,'linewidth',1.1);
 set(gca,'ZGrid','on','GridLineStyle',':','GridColor','k','GridAlpha',0.5);
% line([1 2 ],[xlim],'color','b','LineStyle','-.','linewidth',1.5)
%  
% line([1 2 ],[xlim],'color','b','LineStyle','-.','linewidth',1.5)

 line([0,0],[0,3],[1,1],'color','k','LineStyle','-.','linewidth',1.3)
 line([0,3],[3,3],[1,1],'color','k','LineStyle','-.','linewidth',1.3)
 line([3,3],[3,3],[0,1],'color','k','LineStyle','-.','linewidth',1.3)

 
 
 
