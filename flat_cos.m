%% Flat fading channel, improved Jakes model
% By M. R. Avendi
function h=flat_cos(N,fdTs,tx_pro)
% inputs
% N: numebr of samples
%fdTs: fD*Ts normalized doppler frequency
% R : number of relays in the network, if there is no relay R=0
% output
% h : channel samples

M=18;

% generating uniform random variables
a=-pi;
b=pi;

t=0:N-1;
omega_d=2*pi*fdTs*tx_pro;

Zc=zeros(1,N);
Zs=zeros(1,N);
for n=1:M
    phi_n=a+(b-a).*rand;
    teta=a+(b-a).*rand;
    alfa_n=(2*pi*n-pi+teta)./(4*M);  

    Zc=Zc+sqrt(2/M)*cos(omega_d*t.*cos(alfa_n)+phi_n);
  
    varphi=a+(b-a).*rand;
    Zs=Zs+sqrt(2/M)*cos(omega_d*t.*sin(alfa_n)+varphi);
end    

h=(Zc+1i*Zs)/sqrt(2);

end



