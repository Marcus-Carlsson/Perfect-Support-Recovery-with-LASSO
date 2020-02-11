figure(1);errorbar(linspace(0.07,0.07*max_nl,max_nl),mean(smf_LA_numtrial_nl),std(smf_LA_numtrial_nl),'*b','LineWidth',2);axis([0.05 0.07*max_nl+0.02 -0.5 max(mean(smf_LA_numtrial_nl)+std(smf_LA_numtrial_nl)+2)]);
hold on; errorbar(linspace(0.07+0.005,0.07*max_nl+0.005,max_nl),mean(smf_QEl0_numtrial_nl),std(smf_QEl0_numtrial_nl),'or','LineWidth',2);
errorbar(linspace(0.07-0.005,0.07*max_nl-0.005,max_nl),mean(smf_QEIF_numtrial_nl),std(smf_QEIF_numtrial_nl),'dg','LineWidth',2);
errorbar(linspace(0.07+0.01,0.07*max_nl+0.01,max_nl),mean(smf_Huber_numtrial_nl),std(smf_Huber_numtrial_nl),'xm','LineWidth',2);
errorbar(linspace(0.07-0.01,0.07*max_nl-0.01,max_nl),mean(smf_RWl1_numtrial_nl),std(smf_RWl1_numtrial_nl),'*y','LineWidth',2);
title(['Support misfit for k=',num2str(k)]);
hold off;

figure(2);plot(linspace(0.07,0.07*max_nl,max_nl),mean(dgt_LA_numtrial_nl),'b','LineWidth',2);axis([0.05 0.07*max_nl+0.02 0 0.8]);
hold on;
plot(linspace(0.07,0.07*max_nl,max_nl),mean(dgt_Huber_numtrial_nl),'m','LineWidth',2);
plot(linspace(0.07,0.07*max_nl,max_nl),mean(dgt_RWl1_numtrial_nl),'y','LineWidth',2);
plot(linspace(0.07,0.07*max_nl,max_nl),mean(dgt_QEl0_numtrial_nl),'r','LineWidth',2);
plot(linspace(0.07,0.07*max_nl,max_nl),mean(dgt_QEIF_numtrial_nl),'g','LineWidth',2);
plot(linspace(0.07,0.07*max_nl,max_nl),mean(dgt_OR_numtrial_nl),'k','LineWidth',2);
plot(linspace(0.07,0.07*max_nl,max_nl),mean(dgt_LAPC_numtrial_nl),'b*','LineWidth',2)
title(['Distance to ground truth k=',num2str(k)]);
legend('LASSO','Huber','RWl_1','QEl_0','QEIF','Oracle solution','location','northwest');
hold off;

