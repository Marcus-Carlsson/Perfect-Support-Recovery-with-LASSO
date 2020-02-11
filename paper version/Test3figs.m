%code for generating the figures in test3. The code relies on both Test2
%and Test3 as well as Test2figs have been run and output kept in memory.

%This is the maximum misit between LASSO (x_lambda) and the lambda regularized oracle
%solution (x_lambdaOR)
max(max(max(dlambdaOR_nl_numtrial_r)))

%the below values are mean, std and median for each value of r, using data
%from 10 noiselevels and 20 numtrials, i.e. 200 values. They are used in
%the last plot and also reported in the text

lambda_smf_r=[];
std_lambda_smf_r=[];
median_lambda_smf_r=[];
for k=1:5
    temp=lambda_smf_nl_numtrial_r(:,:,k);
    lambda_smf_r=[lambda_smf_r mean(temp(:))];
    std_lambda_smf_r=[std_lambda_smf_r std(temp(:))];
    median_lambda_smf_r=[median_lambda_smf_r median(temp(:))];
end
lambda_smf_r
std_lambda_smf_r
median_lambda_smf_r

dgt_nl_r=mean(dgt_nl_numtrial_r,2);dgt_nl_r=squeeze(dgt_nl_r);%average over the 20 trials
figure(2);hold on;for k=1:5,plot(dgt_nl_r(:,k),'linewidth',2);end;hold off
figure(3);hold on;for nl=1:10,plot([2:2:10],dgt_nl_r(nl,:)/(0.01*nl),'linewidth',2);end;%the division normalizes by noiselevel
%the below line plots formula (28) using values also from Test 2
plot([2:2:10],(sqrt(1-mean_delta_lower_m_r(5,[1:5])).*(mean_proj_m_r(5,[1:5]))+sqrt(5*[1:5]).*lambda_smf_r./((1-mean_delta_lower_m_r(5,[1:5])))),'*');
hold off



