%This script generates the figures for Test 2, given that Test2.m has been
%run. It focuses on values of sigma_k(A'A), delta_k and the size of the
%projection of a random vector onto the span of a random matrix A with k
%columns, where as usual k=r*m/100

%compute mean and standard deviation for delta_lower
mean_delta_lower_m_r=squeeze(mean(delta_lower_m_r_numtrial,3));
std_delta_lower_m_r=squeeze(std(delta_lower_m_r_numtrial,0,3));

%do the same for the max of delta_lower and delta_upper
delta_m_r_numtrial=max(delta_lower_m_r_numtrial,delta_upper_m_r_numtrial);
mean_delta_m_r=squeeze(mean(delta_m_r_numtrial,3));
std_delta_m_r=squeeze(std(delta_m_r_numtrial,0,3));

%and do it for proj
mean_proj_m_r=squeeze(mean(proj_m_r_numtrial,3));

sizes=[(1:10) (20:10:100) (200:100:400)];
m=50*sizes;

figure(5);hold on;
for k=1:10,
    plot(m,1-mean_delta_lower_m_r(:,k),'linewidth',2);
end
axis([0 max(m) 0 1]);%legend('2','4','6','8','10','12','14','16','18','20');
for k=1:10,
    errorbar(m,1-mean_delta_lower_m_r(:,k),3*std_delta_lower_m_r(:,k))
end
hold off

figure(6);hold on;
for k=1:10,
    plot(m,mean_delta_m_r(:,k),'linewidth',2);
end
axis([0 max(m) 0 1.2]);%legend('1','2','3','4','5','6','7','8','9','10');
for k=1:10,
    errorbar(m,mean_delta_m_r(:,k),3*std_delta_m_r(:,k))
end
hold off

%The last graph was not used..., but the values are still used in Test3figs
%to check validity of formula (25). To access the values for m=250 (reported in the paper), type
%mean_proj_m_r(9,:)

figure(7);hold on;
for k=1:10,
    plot(m,mean_proj_m_r(:,k),'linewidth',2);
end
hold off