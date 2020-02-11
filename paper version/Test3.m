%This script was used to generate the data in Test3 of the paper
%Perfect support recovery with LASSO. The script generats simulations over
%3 different parameters, noiselevel (abbreviated nl), r (ratio k/m in percent) and numtrial which
%simply repeats the experiment 20 times
%The parameters m (amount of rows) and n (amount of columns) are held fixed



%the following are the output arrays
lambda_smf_nl_numtrial_r=[];%stores optimal value of lambda for support of x_lambda to be inside support xgt
dgt_nl_numtrial_r=[];%distance to ground truth depending on nl (noiselevel) numtrial and r (for the lambda found above)
dlambdaOR_nl_numtrial_r=[];%distance to lambdaregularized oracle solution (for the lambda found above)

minnoiselevel=0.01;%set to 0.01 to run noise between 1 and 10%
m=250;%number of rows
n=2*m;%twice as many columns as rows


for noisecount=1:10,%noiselevels between minnoiselevel and 10*minnoiselevel
    nl=noisecount*minnoiselevel
    dgt_numtrial_r=[];%These are internal variables which will be used to construct the main ones with same name above
    dlambdaOR_numtrial_r=[];
    lambda_smf_numtrial_r=[];
    
    
for numtrial=1:20 %we run 20 experiments for each nl and each r (r loop comes later)
A=randn(m,n);e=randn(m,1);%creating real data, next line complex, so comment next one out if you wanna run real data
A=randn(m,n)+i*randn(m,n);e=randn(m,1)+i*randn(m,1);
A=A./vecnorm(A);e=nl*e/norm(e);%normalize columns of A, set noise level

AtA=A'*A;%for faster evaluation we store this
t=0.9/norm(A)^2;%stepsize parameter in FBS

dgt_r=[];
dlambdaOR_r=[];
lambda_smf_r=[];

for k=5:5:25 %here comes the r loop, running values of k from 2 percent to 10 percent (of m)
xgt=zeros(n,1);   %we construct ground truth, again comment out second row if you want to work with real data
xgt(1:k,1)=sign(randn).*(1+rand(k,1));%max ratio 2 of max and min nonzero value in xgt
xgt=zeros(n,1);xgt(1:k,1)=exp(2*pi*1i*rand(k,1)).*(1+rand(k,1));
xgt=xgt/norm(xgt);%normalize it

b=A*xgt+e;%measurement
Atilde=A(:,[1:k]);xo=zeros(n,1);xo(1:k)=lsqminnorm(Atilde,b); %create oracle solution


dgt=1;%set dgt to 1 
dlambdaOR=1; %if the value stays one then correct support was never found
lambda_smf=0;

lambdasteps=[0.03*[1:25] 0.15*[6:30] 0.3*[16:20]];

%we compute LASSO solution using FBS on the dual formulation for various
%values of lambda (i.e. we minimize lambda||x||_1+1/2||Ax-b||^2)
for lambda=nl*lambdasteps%this fine grid of lambdas start with 5% of noiselevel, which is too little for correct 
    %support retrieval at least for noise above 1%, and takes 5% steps up
    %to 1.25*noiselevel. This is almost always sufficient to find correct
    %support, but sometimes even for small noise (<10%) and small k (<5% of
    %n), higher values up to 10*nl is needed. For this reason, the stepsize
    %increases by 50% steps up to 15*nl
    %The scale is here more coarse than in previous tests becuase much
    %higher values of lambda are needed for correct support recovery than
    %just minimizing dgt

    %we first compute the lambda regularized oracle solution
    x_lambdaORtrunc=zeros(k,1);%initialize x lambda oracle, where we only operate over the active first k
    for j=1:2000
        x_lambdaORtrunc=proxL1(x_lambdaORtrunc-t*Atilde'*(Atilde*x_lambdaORtrunc-b),t*lambda);
    end
    x_lambdaOR=zeros(n,1);x_lambdaOR(1:k)=x_lambdaORtrunc;

    
    convergedLASSO=0;%this parameter is commented out, uncomment if u wanna check that the routine converges before max number of j
    x_lambda=zeros(n,1);%initialize xLASSO
    for j=1:10000
        x_lambda_previous=x_lambda;%store previous value
        x_lambda=proxL1(x_lambda-t*(AtA*x_lambda-A'*b),t*lambda);%this is the FBS update step
        if norm(x_lambda-x_lambda_previous)<10^(-16)%break the loop if increments are too small, i.e. we consider this as convergence
            break
        else
            convergedLASSO=convergedLASSO+1;%counts the number of steps until algorithm converges
        end
    end %uncomment below section if you want to see what is happening
%    lambda
%    convergedLASSO
%    figure(1);
%    plot(log10(abs(xlambdaOR)+10^(-10)),'r');%axis([1 n 0 1.2*norm(xo,inf)]); 
%    title('LASSO'); hold on; plot(log10(abs(xLASSO)+10^(-10)));legend({'lambda regularized oracle solution','reconstruction for different \lambda'}); hold off 
 %   pause(1.3)

 %security check that the algorithm converged
     if convergedLASSO == 10000,
         disp('Not enough iterations')
         convergedLASSO
         lambda/noiselevel
     end
 
    %we decide that correct support has been found once the "off support
    %norm" of x_lambda is small enough. In our experience, it makes no
    %difference with checking when it is identically zero
    if norm(x_lambda(k+1:n)) < 10^(-8),
        dgt=norm(xgt-x_lambda);%store distance to ground truth for smallest value of lambda that retrieves correct support (or a subset thereof)
        lambda_smf=lambda/nl;%as well as the lambda for which it happened
        dlambdaOR=norm(x_lambdaOR-x_lambda);
          break;
    end
end %lambda loop ends
dgt_r=[dgt_r dgt];
lambda_smf_r=[lambda_smf_r lambda_smf];
dlambdaOR_r=[dlambdaOR_r dlambdaOR];

end %kloop ends
dgt_numtrial_r=[dgt_numtrial_r;dgt_r];
lambda_smf_numtrial_r=[lambda_smf_numtrial_r;lambda_smf_r];
dlambdaOR_numtrial_r=[dlambdaOR_numtrial_r;dlambdaOR_r];

end %numtrial loop ends

dgt_nl_numtrial_r(noisecount,:,:)=dgt_numtrial_r;
lambda_smf_nl_numtrial_r(noisecount,:,:)=lambda_smf_numtrial_r;
dlambdaOR_nl_numtrial_r(noisecount,:,:)=dlambdaOR_numtrial_r;
end %loop over nl ends

%This ends the simulation, left to do is the related figures. The script
%for this is found in separate file called Test2figs