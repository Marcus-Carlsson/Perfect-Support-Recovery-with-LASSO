%This script was used to generate the data in Test 1 of the paper
%Perfect support recovery with LASSO. The script generats simulations over
%3 different parameters, m (amount of rows), r (relative amount of non-zeros in
%ground truth) and numtrial (which simply repeats the experiment for same
%value of m and r). Noiselevel is fixed at 0.1%

%Activate plot on row 82 to see what is going on...

noiselevel=0.1;

max_m=19;%Should be 19 to be compatible with Test1figs 
max_numtrial=20;
max_r=10;%Should be 10 to be compatible with Test1figs 

%the following are the output arrays
r_m=zeros(1,max_m);%the maximum value of r (i.e. 100*k/n) which can be retrieved in 100% of the trials
smf_m_numtrial_r=zeros(max_m,max_numtrial,max_r);%support misfit as multiarray depending on m, numtrial and r, in that order
dgt_m_numtrial_r=zeros(max_m,max_numtrial,max_r);%distance to ground truth 
dOR_m_numtrial_r=zeros(max_m,max_numtrial,max_r);%distance to oracle solution (not used in the plots)
M_m_numtrial=zeros(max_m,max_numtrial,max_r);%records mutual coherence number as function of n and numtrial (not used in plots)
lambda_smf_m_numtrial_r=zeros(max_m,max_numtrial,max_r);%stores optimal value of lambda for best smf, divided by noiselevel
lambda_dgt_m_numtrial_r=zeros(max_m,max_numtrial,max_r);%stores optimal value of lambda for dgt, divided by noiselevel


parfor numtrial=1:max_numtrial %we run 20 experiments for each n and each r (r loop comes later)

numberofrowsfirstexperiment=25;%should be divisible by 25
for iter_m=1:max_m,
if iter_m<10 %lets m run from 25 to 225 with increments of 25 and then from 250 to 2500 with increments of 250
    m=numberofrowsfirstexperiment*iter_m;
else
    m=numberofrowsfirstexperiment*10*(iter_m-9);
end
    n=2*m%twice as many columns as rows
    
A=randn(m,n);e=randn(m,1);%creating real data, next line complex, so comment next one out if you wanna run real data
A=randn(m,n)+i*randn(m,n);e=randn(m,1)+i*randn(m,1);
A=A./vecnorm(A);e=noiselevel*e/norm(e);%normalize columns of A, set noise level
M_m_numtrial(iter_m,numtrial)=max(max(abs(A'*A)-eye(n,n)));%computes the mutual coherence number of A

t=0.9/norm(A)^2;%stepsize parameter in FBS (theoretical upper bound is 1/norm(A)^2)

for iter_r=1:max_r %here comes the r loop, running values of k for r from 4 percent to 40 percent 
    k=m/25*iter_r;
xgt=zeros(n,1);   %we construct ground truth, again comment out second row if you want to work with real data
xgt(1:k,1)=sign(randn).*(1+rand(k,1));%max ratio 2 of max and min nonzero value in xgt
xgt=zeros(n,1);xgt(1:k,1)=exp(2*pi*1i*rand(k,1)).*(1+rand(k,1));
xgt=xgt/norm(xgt);%normalize it



b=A*xgt+e;%create measurement
Atilde=A(:,[1:k]);xo=zeros(n,1);xo(1:k)=lsqminnorm(Atilde,b); %create oracle solution

AtA=A'*A;%these are stored for speed
Atb=A'*b;

smf=n;%set smf to maxvalue
dgt=norm(xgt);%set dgt to 1 
lambda_smf=0;
lambda_dgt=0;

lambdasteps=[0.01*[1:25] 0.05*[6:30] 0.1*[16:50]];

%we compute LASSO solution using FBS on the dual formulation for various
%values of lambda (i.e. we minimize lambda||x||_1+1/2||Ax-b||^2)
for lambda=noiselevel*lambdasteps%this fine grid of lambdas start with 1% of noiselevel, which is too little for correct 
    %support retrieval at least for noise above 1%, and takes 1% steps up
    %to 0.25*noiselevel. This is almost always sufficient to find correct
    %support, but sometimes even for small noise (<10%) and small r (<10% of
    %m), higher values are needed. For this reason, the stepsize
    %increases by 10% steps up to 5*nl, by which time the smf function is
    %always increasing so no need to try larger lambdas
    convergedLASSO=0;%this parameter is commented out, uncomment if u wanna check that the routine converges before max number of j
    x_lambda=zeros(n,1);%initialize x_lambda
    x_lambda_previous=x_lambda;
    for j=1:100000
        %due to the extremely slow convergence of FBS for small lambdas, we first do a FISTA
        %iteration, but switch to FBS either after 500 iterations or when
        %the iterates get sufficiently close. We do not use FBS in the end
        %due to its oscillatory behavior. The delimiters 500 and 10^(-4) are a chosen a bit ad hoc, but makes sure that in the end the 
        %FBS and not FISTA is in use.
        if convergedLASSO<500 & norm(x_lambda-x_lambda_previous)>10^(-4),
        x_temp=x_lambda+(j-1)/(j+2)*(x_lambda-x_lambda_previous);
        x_lambda_previous=x_lambda;%store previous value
        x_lambda=proxL1(x_temp-t*(AtA*x_temp-Atb),t*lambda);
        %this is the standard FBS update step
        else
            x_lambda_previous=x_lambda;
            x_lambda=proxL1(x_lambda-t*(AtA*x_lambda-Atb),t*lambda);%this is the FBS update step
               end
        if norm(x_lambda-x_lambda_previous)<10^(-8)%break the loop if increments are too small, i.e. we consider this as convergence
            %For lambda around 0.01*noiselevel the code is very sloow, and
            %may need 200000 iterations to converge to 10^(-16) (whereas above 0.1 it is usually a vew 100 iterations only). In a
            %separate script we verified that all output parameters are the
            %same if we stop at 10^(-8), and therefore we chose this as stopping criteria which
            %reduces computational times
            break
        else
            convergedLASSO=convergedLASSO+1;%counts the number of steps until algorithm converges
        end
    end
%    convergedLASSO
    %This is a security warning that fires only if the algorithm did not
    %terminate before 100000 iterations
    if convergedLASSO == 100000,
         disp('Not enough iterations')
         convergedLASSO
         lambda/noiselevel
    end
    
%    Uncomment this block if you wanna see what is going on
%   lambda
%   convergedLASSO
%   figure(1);%this allows you to monitor what is going on, comment out for faster evaluation
%   plot(abs(xo),'r');axis([1 n 0 1.2*norm(xo,inf)]); title('LASSO'); hold on; plot(abs(x_lambda));legend({'oracle solution','reconstruction for different \lambda'}); hold off 

    supp=zeros(n,1);id=find(abs(x_lambda)>10^(-8));supp(id)=1;%create vector with ones on the support of x_lambda,
    %here we could just as well put ~=0, but we do >10^(-8) to avoid
    %the possibility of missing zeroes due to numerical errors. In our
    %experience, there is no difference
    truesupp=ones(k,1);truesupp=[truesupp;zeros(n-k,1)];%create vector with ones on the support of xgt
    if norm(truesupp-supp,1)<smf,
        smf=min(smf,norm(truesupp-supp,1));%store best value of support misfit
        lambda_smf=lambda/noiselevel;%as well as the lambda for which it happened, divided by noiselevel
    end
    if norm(xgt-x_lambda)<dgt,
        dgt=min(dgt,norm(xgt-x_lambda));%store best value of distance to ground truth
        lambda_dgt=lambda/noiselevel;%as well as the lambda for which it happened, divided by noiselevel
    end
%if lambda is large enough that support of x_lambda is subset of xgt, no need to continue increasing lambda
    if all(ismember(id,[1:k])) 
        break, 
    end
end %lambda loop ends

%sometimes 1% of the noise still gives you the best fit for dgt. However,
%going lower makes the algorithm so slow that it is not feasible to
%continue shrinking lambda. The limit case lambda=0 is the hard constrained
%l1-problem, also known as basis pursuit. We therefore compute this
%solution separately with an ADMM solver.

rho=100;
convergedLASSO2=0;
x_0=zeros(n,1);
z=zeros(n,1);
u=zeros(n,1);
AAt = A*A';
P = eye(n) - A' * (AAt \ A);
q = A' * (AAt \ b);
        x_0old=ones(n,1);
    for j=1:100000
    x_0=proxL1(z+u,1/rho);
    z=P*(x_0-u)+q;
    u=u+z-x_0;
        if norm(x_0-x_0old)<10^-(8)
            break
        else
            convergedLASSO2=convergedLASSO2+1;
        end
                x_0old=x_0;
    end
     if convergedLASSO2 == 100000,
         disp('Not enough iterations basis pursuit')
     end
     supp=zeros(n,1);id=find(abs(x_0)>10^(-8));supp(id)=1;
if norm(truesupp-supp,1)<smf,%norm(truesupp-supp,1) computes the amount of wrong non-zeroes and wrong zeroes
    smf=norm(truesupp-supp,1);
    lambda_smf=0;
end
if norm(xgt-x_0)<dgt,
    dgt=norm(xgt-x_0);
    lambda_dgt=0;
end







smf_m_numtrial_r(iter_m,numtrial,iter_r)=smf;%we store the optimal values for this value of k
lambda_smf_m_numtrial_r(iter_m,numtrial,iter_r)=lambda_smf;
dgt_m_numtrial_r(iter_m,numtrial,iter_r)=dgt;
lambda_dgt_m_numtrial_r(iter_m,numtrial,iter_r)=lambda_dgt;
dOR_m_numtrial_r(iter_m,numtrial,iter_r)=norm(xgt-xo);

end %kloop ends
end %loop over m ends

end %numtrial loop ends

%this final loop computes the values for perfect support recovery graph,
%figure 3

for iter_m=1:max_m,
id=find(vecnorm(squeeze(smf_m_numtrial_r(iter_m,:,:)))==0);%stores indices such that smf was 0 for all numtrials
if size(id) == 0,rmax=0;
else rmax=max(id);end%let rmax be the maximum such index

r_m(iter_m)=rmax;%add this value of r to r_m
end

%This ends the simulation, left to do is the related figures. The script
%for this is found in separate file called Test1figs