%code for generating the figures in test1, given that r_m, smf_m_numtrial_r and
%dgt_m_numtrial_r have been created with Test1.m

%first figure (K plot)
m=[[1:9]'*25;[1:10]'*250];n=2*m;%create a vector containing sizes of m and n used
figure(1);plot(m,4*r_m,'linewidth',2);%a value r=1 corresponds to 4 percent
title('Perfect support recovery');
axis([0 2500 0 13]);

%second figure (smf plot)
smf_m_r=mean(smf_m_numtrial_r,2);smf_m_r=squeeze(smf_m_r);%average over the 20 trials
figure(2);hold on;for k=1:10,plot(m,[100*smf_m_r(:,k)./(4*k*m/100)],'linewidth',2);end;hold off
title('Support misfit as percentage of n');

%third figure (dgt plot)
dgt_m_r=mean(dgt_m_numtrial_r,2);dgt_m_r=squeeze(dgt_m_r);
figure(3);hold on;plot([4:4:40],[dgt_m_r(1,:)],'linewidth',2);
plot([4:4:40],[dgt_m_r(10,:)],'linewidth',2);
plot([4:4:40],[dgt_m_r(19,:)],'linewidth',2);
hold off
title('Distance to ground truth');
legend('m=25','m=250','m=2500');

%alternative full (dgt plot)
dgt_m_r=mean(dgt_m_numtrial_r,2);dgt_m_r=squeeze(dgt_m_r);
figure(4);hold on;for k=1:10,
    plot([4:4:40],[dgt_m_r(k,:)],'linewidth',2);
end
hold off
title('Distance to ground truth');


