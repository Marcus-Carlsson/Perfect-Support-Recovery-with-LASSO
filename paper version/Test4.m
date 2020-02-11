%This script launches simulations with the compressed sensing standard
%setup Ax=b with a mxn matrix A and normalized columns, vector x of cardinality
%k and Gaussian noise of magnitude epsilon, which was used in test 4 of the
%paper "perfect support recovery with LASSO" by M. Carlsson
%We compare the following routines:
%dual LASSO
%dual LASSO with post correction
%Huber penalty by I. Selesnick
%Reweighted least squares by Candes, Wakin, Boyd
%Quadratic Envelope of mu times cardinality as well as the indicator function of k non-zeros, by Carlsson, Olsson, Gerosa

m=80;%matrix rows
n=200;%matrix columns
k=15;%amount of non-zeroes

levelstoppingcriteria=10^(-8);%accuracy of two consequtive iterates for stopping FBS
maxiter=100000;%maximum iterations if above stopping does not happen, this will flash a warning

%The script runs various tests for various noiselevels,
numtrial=20;%number of trials for each noise level
max_nl=5;%number of noise levels

smf_LA_numtrial_nl=zeros(numtrial,max_nl);%best value of support misfit for LASSO
smf_Huber_numtrial_nl=zeros(numtrial,max_nl);%same for the Huber method
smf_RWl1_numtrial_nl=zeros(numtrial,max_nl);%same for the reweighted least squares method
smf_QEl0_numtrial_nl=zeros(numtrial,max_nl);%same for QEl0
smf_QEIF_numtrial_nl=zeros(numtrial,max_nl);%same for QEIF

%the next block creates similar variables for distance to ground truth,
%with the addition of dgt_LAPCtot which is the distance to ground truth
%from the LASSO post corrected method, (i.e. solving a least squares
%problem using the columns of A found by the support of the LASSO solution), 
%as well as the distance to ground truth from the oracle 
dgt_LA_numtrial_nl=zeros(numtrial,max_nl);
dgt_Huber_numtrial_nl=zeros(numtrial,max_nl);
dgt_RWl1_numtrial_nl=zeros(numtrial,max_nl);
dgt_QEl0_numtrial_nl=zeros(numtrial,max_nl);
dgt_QEIF_numtrial_nl=zeros(numtrial,max_nl);
dgt_LAPC_numtrial_nl=zeros(numtrial,max_nl);
dgt_OR_numtrial_nl=zeros(numtrial,max_nl);

%the following 3 blocks are just control variables, to test for which
%lambda (or equivalent parameter) divided by noise level, the best smf or
%dgt happened, and how many iterations it took. These are not used in the
%plots but serve to check that we have found a good range of parameter
%values, and can be used to test good parameter regimes...

lambda_smf_LA_numtrial_nl=zeros(numtrial,max_nl);%this stores the lambda value (divided by noise) for which best smf happened
lambda_smf_Huber_numtrial_nl=zeros(numtrial,max_nl);%....
lambda_smf_RWl1_numtrial_nl=zeros(numtrial,max_nl);%....
sqrtmu_smf_QEl0_numtrial_nl=zeros(numtrial,max_nl);%
gamma_smf_QEIF_numtrial_nl=zeros(numtrial,max_nl);%

lambda_dgt_LA_numtrial_nl=zeros(numtrial,max_nl);%does the same for dgt
lambda_dgt_Huber_numtrial_nl=zeros(numtrial,max_nl);%
lambda_dgt_RWl1_numtrial_nl=zeros(numtrial,max_nl);%
sqrtmu_dgt_QEl0_numtrial_nl=zeros(numtrial,max_nl);%
gamma_dgt_QEIF_numtrial_nl=zeros(numtrial,max_nl);%


%Since the two quadratic envelope methods are capable of finding the oracle
%solution, these two additional values register for which parameters this
%happens

sqrtmu_foundoracle_QEl0tot=[];%this stores sqrtmu/nl
gamma_foundoracle_QEIFtot=[];%this simply stores gamma

for l=1:numtrial,
tic
%here noiselevel steps are set at 0.07, which of course can be modified at
%wish, but it requires manual update of the file Test4fig.m
for iter_nl=1:max_nl,
    
    nl=0.07*iter_nl;
    
A=randn(m,n);e=randn(m,1);xgt=zeros(n,1);xgt(1:k,1)=sign(randn).*(1+rand(k,1));xgt=xgt/norm(xgt);%creating real data, next line complex, so comment next one out if you wanna run real data
A=randn(m,n)+i*randn(m,n);e=randn(m,1)+i*randn(m,1);xgt=zeros(n,1);xgt(1:k,1)=sign(randn).*(1+rand(k,1))+i*sign(randn).*(1+rand(k,1));xgt=xgt/norm(xgt);%create ground truth
A=A./vecnorm(A);e=nl*e/norm(e);%normalize columns of A, set noise level
b=A*xgt+e;%measurement
Atilde=A(:,[1:k]);xo=zeros(n,1);xo(1:k)=lsqminnorm(Atilde,b); %create oracle solution
t=0.9/norm(A)^2;%stepsize in FBS
AtA=A'*A;%stored for computational speed
truesupp=ones(k,1);truesupp=[truesupp;zeros(n-k,1)];%true support of ground truth

dgt_OR_numtrial_nl(l,iter_nl)=norm(xgt-xo);%first output variable, distance from ground truth to oracle solution

%we compute LASSO solution using FBS on the dual formulation for various
%values of lambda (i.e. we minimize lambda||x||_1+1/2||Ax-b||^2)

%temporary variables to store best values:
smf_LA=n;%best support misfit
dgt_LA=1000;%best distance to ground truth, 1000 is just an adhoc value large enough to be bigger than any output
dgt_LAPC=1000;%best distance to groung truth with post correction
maxconverged_LA=0;

for lambda=nl*[0.01*[1:25] 0.05*[6:30] 0.1*[16:30]]%By testing we have found this range of values good. 
%    If best values occur for minimumor maximum lambda a warning will go
%    off
    supp=zeros(n,1);
    convergedLASSO=0;
    x_lambda=zeros(n,1);%initialize x_lambda
    for j=1:maxiter
        x_lambda_previous=x_lambda;%store previous value
        x_lambda=proxL1(x_lambda-t*(AtA*x_lambda-A'*b),t*lambda);
        if norm(x_lambda-x_lambda_previous)<levelstoppingcriteria
            break
        else
            convergedLASSO=convergedLASSO+1;
        end
    end
     if convergedLASSO == maxiter,
         disp('Not enough iterations LASSO')
     end
%   figure(1);
%   plot(abs(xo),'r');axis([1 n 0 1.2*norm(xo,inf)]); title('LASSO'); hold on; plot(abs(x_lambda));legend({'oracle solution','reconstruction for different \lambda'}); hold off 
%   pause(0.3)
id=find(abs(x_lambda)>10^(-7));supp(id)=1;%here we could just as well put ~=0, but we do >10^(-7) to avoid
%the possibility of missing zeroes due to numerical errors

%store the best value for the best lambda/nl (not necessarily same for smf and dgt)

if norm(truesupp-supp,1)<smf_LA,%norm(truesupp-supp,1) computes the amount of wrong non-zeroes and wrong zeroes
    smf_LA=norm(truesupp-supp,1);
    lambda_smf_LA=lambda/nl;
end
if norm(xgt-x_lambda)<dgt_LA,
    dgt_LA=norm(xgt-x_lambda);
    lambda_dgt_LA=lambda/nl;
end
Atemp=A(:,id);xtemp=zeros(n,1);xtemp(id)=lsqminnorm(Atemp,b); %post correct using support of x_lambda
dgt_LAPC=min(dgt_LAPC,norm(xgt-xtemp));


%if lambda is large enough that support of x_lambda is subset of xgt, or the current iterate is worse than 2 times best,
%no need to continue increasing lambda
    if all(ismember(id,[1:k])) 
        break, 
    end
end

%warnings for poor parameter ranges
if lambda_smf_LA == 3
    disp('Warning, poor parameter range for LASSO smf')
    lambda_smf_LA
    iter_nl
end
if lambda_dgt_LA == 3
    disp('Warning, poor parameter range for LASSO dgt')
    lambda_dgt_LA
    iter_nl
end


%finally we solve the problem with lambda=0, i.e. the hard constrained one,
%using instead ADMM, (which converges 100 times faster at least than FBS with lambdas of size 1% or less)
rho=100;
convergedLASSO2=0;
x_0=zeros(n,1);
z=zeros(n,1);
u=zeros(n,1);
AAt = A*A';
P = eye(n) - A' * (AAt \ A);
q = A' * (AAt \ b);
        x_0old=ones(n,1);
    for j=1:maxiter
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
     if convergedLASSO2 == maxiter,
         disp('Not enough iterations basis pursuit')
     end
     id=find(abs(x_0)>10^(-7));supp(id)=1;
if norm(truesupp-supp,1)<smf_LA,%norm(truesupp-supp,1) computes the amount of wrong non-zeroes and wrong zeroes
    smf_LA=norm(truesupp-supp,1);
    lambda_smf_LA=0;
end
if norm(xgt-x_0)<dgt_LA,
    dgt_LA=norm(xgt-x_0);
    lambda_dgt_LA=0;
end
Atemp=A(:,id);xtemp=zeros(n,1);xtemp(id)=lsqminnorm(Atemp,b); %post correct using support of x_lambda
dgt_LAPC=min(dgt_LAPC,norm(xgt-xtemp));




%we finally store best values for this noiselevel and numtrial
smf_LA_numtrial_nl(l,iter_nl)=smf_LA;
dgt_LA_numtrial_nl(l,iter_nl)=dgt_LA;
dgt_LAPC_numtrial_nl(l,iter_nl)=dgt_LAPC;
lambda_smf_LA_numtrial_nl(l,iter_nl)=lambda_smf_LA;
lambda_dgt_LA_numtrial_nl(l,iter_nl)=lambda_dgt_LA;




%We now run the algorithm for Huber, following precisely the setup in
%Selesnicks paper sparse regularization via convex analysis. The structure
%of the code is a copy paste of the above

smf_Huber=n;
dgt_Huber=1000;

for lambda=nl*[0.03*[1:50] 0.1*[16:30]]%By testing we have found this range of values good. 
    %If best values occur for minimumor maximum lambda a warning will go
    %off
    convergedHuber=0;
    supp=zeros(n,1);
    xHuber=zeros(n,1);%initialize xHuber
    vHuber=zeros(n,1);
    uHuber=zeros(n,1);
    wHuber=zeros(n,1);
    gamma=0.8;%see Selesnicks paper for this parameter choise
    rho=gamma/(1-gamma)*norm(AtA);
    mu=1/rho;
    for j=1:maxiter
        xHuberPrevious=xHuber;
        wHuber=xHuber-mu*(AtA*(xHuber+gamma*(vHuber-xHuber))-A'*b);
        uHuber=vHuber-mu*gamma*AtA*(vHuber-xHuber);
        xHuber=proxL1(wHuber,mu*lambda);
        vHuber=proxL1(uHuber,mu*lambda);
        if norm(xHuber-xHuberPrevious)<levelstoppingcriteria
            break
        else
            convergedHuber=convergedHuber+1;
        end
    end
         if convergedHuber == maxiter
         disp('Not enough iterations Huber')
     end

%     lambda
%     convergedHuber
%    figure(2);
%    plot(abs(xo),'r');axis([1 n 0 1.2*norm(xo,inf)]); title('Huber'); hold on; plot(abs(xHuber));legend({'oracle solution','reconstruction for different \lambda'}); hold off 
%    pause(0.3)
id=find(abs(xHuber)>10^(-7));supp(id)=1;
if norm(truesupp-supp,1)<smf_Huber,
    smf_Huber=norm(truesupp-supp,1);
    lambda_smf_Huber=lambda/nl;
end
if norm(xgt-xHuber)<dgt_Huber,
    dgt_Huber=norm(xgt-xHuber);
    lambda_dgt_Huber=lambda/nl;
end
    if all(ismember(id,[1:k]))
        break, 
    end%if lambda is large enough that support of xHuber is subset of xgt, no need to continue increasing lambda
end

if lambda_smf_Huber == 0.03 || lambda_smf_Huber == 3
    disp('Warning, poor parameter range for Huber smf')
    lambda_smf_Huber
    iter_nl
end
if lambda_dgt_Huber == 0.03 || lambda_dgt_Huber == 3
    disp('Warning, poor parameter range for Huber dgt')
    lambda_dgt_Huber
    iter_nl
end


smf_Huber_numtrial_nl(l,iter_nl)=smf_Huber;
dgt_Huber_numtrial_nl(l,iter_nl)=dgt_Huber;
lambda_smf_Huber_numtrial_nl(l,iter_nl)=lambda_smf_Huber;
lambda_dgt_Huber_numtrial_nl(l,iter_nl)=lambda_dgt_Huber;






%We now run the algorithm Reweighted Least Squares, following the paper by
%Candes, Wakin, Boyd as closely as possible

smf_RWl1=n;
dgt_RWl1=1000;
maxconverged_RWl1=0;
epsilon=0.0001;%this parameterchoise is recommended somewhere in the litterature....
xRWl1old=ones(n,1);
for lambda=nl*[0.01*[3:30]]% this algorithm is more sensitive to lambda, and dont need as high values of lambda.
    supp=zeros(n,1);
    xRWl1=zeros(n,1);%initialize xRWLS
    w=lambda*ones(n,1);%set weights to 1
    for updateweights=1:10%this is the extra loop for this algorithm. 10 is chosen ad hoc, but by visual 
        %inspection it seems it is enough for the algorithm to converge
    convergedRWl1=0;
       for j=1:maxiter
           xRWl1_previous=xRWl1;
           xRWl1=proxL1(xRWl1-t*(AtA*xRWl1-A'*b),t*w);
           if norm(xRWl1-xRWl1_previous)<levelstoppingcriteria
               break
           else
               convergedRWl1=convergedRWl1+1;
           end
        end
        if convergedRWl1 == maxiter
                  disp('Not enough iterations RWLS')
     end
     w=lambda^2./(abs(xRWl1)+epsilon);
 %  figure(2);
 %  plot(abs(xo),'r');axis([1 n 0 1.2*norm(xo,inf)]); title('RWl1'); hold on; plot(abs(xRWl1));legend({'oracle solution','reconstruction for different \lambda'}); hold off 
 %   pause(0.3)
   if norm(xRWl1-xRWl1old)<levelstoppingcriteria;%break outer loop is no changes happened
       break;
   end     
   xRWl1old=xRWl1;
end
id=find(abs(xRWl1)>10^(-7));supp(id)=1;
if norm(truesupp-supp,1)<smf_RWl1,
    smf_RWl1=norm(truesupp-supp,1);
    lambda_smf_RWl1=lambda/nl;
end
if norm(xgt-xRWl1)<dgt_RWl1,
    dgt_RWl1=norm(xgt-xRWl1);
    lambda_dgt_RWl1=lambda/nl;
end
    if all(ismember(id,[1:k]))
        break, 
    end%if lambda is large enough that support of xRWLS is subset of xgt, no need to continue increasing lambda
end
if lambda_smf_RWl1 == 0.03 || lambda_smf_RWl1 == 0.3
    disp('Warning, poor parameter range for RWl1 smf')
    lambda_smf_RWl1
    iter_nl
end
if lambda_dgt_RWl1 == 0.03 || lambda_dgt_RWl1 == 0.3
    disp('Warning, poor parameter range for RWl1 dgt')    
    lambda_dgt_RWl1
    iter_nl
end


smf_RWl1_numtrial_nl(l,iter_nl)=smf_RWl1;
dgt_RWl1_numtrial_nl(l,iter_nl)=dgt_RWl1;
lambda_smf_RWl1_numtrial_nl(l,iter_nl)=lambda_smf_RWl1;
lambda_dgt_RWl1_numtrial_nl(l,iter_nl)=lambda_dgt_RWl1;




%Now we run the method QEl0 from Carlsson, Geroasa, Olssons paper An
%unbiased approach to compressed sensing, with the difference that the
%gamma parameter is set to 0.7 rather than 1 as in that paper (simply because it seems to work better by empirical evidence)
%One could of course run a loop over different values of gamma as well, but
%since anyways this outperforms the above methods most of the time, we keep
%it simple

smf_QEl0=n;
dgt_QEl0=1000;
gamma=0.7;
%these two are new variables which register when the oracle solution are
%found
    sqrtmu_foundoracle_QEl0=zeros(1,40);
    count=0;
%for gamma=0.3:0.2:0.7
    for sqrtmu=nl*[0.03:0.03:0.39] %this grid is very coarse compared with the previous, 
    %but the algorithm usually finds the correct point anyways.
    %We run over parameter sqrtmu, square root of mu, since this has a
    %similar role in the proximal operator as lambda in soft thresholding.
    supp=zeros(n,1);
    convergedQEl0=0;
    x_QEl0=zeros(n,1);%initialize xQEl0
    for j=1:maxiter
        x_QEl0_previous=x_QEl0;
        x_QEl0=ProxQmucard(x_QEl0-t*(AtA*x_QEl0-A'*b),sqrtmu^2,gamma,1/t);
        if norm(x_QEl0-x_QEl0_previous)<levelstoppingcriteria
            break
        else
            convergedQEl0=convergedQEl0+1;
        end
    end
    if convergedQEl0 == maxiter,
         disp('Not enough iterations QEl0')
    end
    id=find(abs(x_QEl0)>10^(-7));supp(id)=1;
%     figure(2);
%     plot(abs(xo),'r');axis([1 n 0 1.2*norm(xo,inf)]);title('QEl0'); hold on; plot(abs(x_QEl0));legend({'oracle solution','reconstruction for different \mu'}); hold off 
%     pause(0.3)

%this one is new compared to the previous 3, and records for which sqrtmu/nl the method
%finds the oracle solution
if norm(xo-x_QEl0)<100*levelstoppingcriteria,
    count=count+1;
    sqrtmu_foundoracle_QEl0(count)=sqrtmu/nl;
end

if norm(truesupp-supp,1)<smf_QEl0,
    smf_QEl0=norm(truesupp-supp,1);
    sqrtmu_smf_QEl0=sqrtmu/nl;
end
if norm(xgt-x_QEl0)<dgt_QEl0,
    dgt_QEl0=norm(xgt-x_QEl0);
    sqrtmu_dgt_QEl0=sqrtmu/nl;
end
%we dont have a stopping criterion here since anyways we test over very few
%values of sqrtmu
end
%end
if (sqrtmu_smf_QEl0 == 0.03 || sqrtmu_smf_QEl0 == 0.039) & count == 0,
    disp('Warning, poor parameter range for QEl0 smf')
    sqrtmu_smf_QEl0
    iter_nl
end
if (sqrtmu_dgt_QEl0 == 0.03 || sqrtmu_dgt_QEl0 == 0.039) & count == 0,
    disp('Warning, poor parameter range for QEl0 dgt')
    sqrtmu_dgt_QEl0
    iter_nl
end

smf_QEl0_numtrial_nl(l,iter_nl)=smf_QEl0;
dgt_QEl0_numtrial_nl(l,iter_nl)=dgt_QEl0;
sqrtmu_smf_QEl0_numtrial_nl(l,iter_nl)=sqrtmu_smf_QEl0;
sqrtmu_dgt_QEl0_numtrial_nl(l,iter_nl)=sqrtmu_dgt_QEl0;
sqrtmu_foundoracle_QEl0tot=[sqrtmu_foundoracle_QEl0tot;sqrtmu_foundoracle_QEl0];



%Finally we run QEIF. Since this uses information about k, it is much
%simpler and the only parameter is gamma, which we run for ~10 different
%values between 0.1 and 1.


smf_QEIF=n;
dgt_QEIF=1000;
     gamma_foundoracle_QEIF=zeros(1,12);
     count=0;

for gamma=1:-0.1:0.1
    supp=zeros(n,1);
    convergedQEIF=0;
    xQEIF=zeros(n,1);%initialize xQEIF
    for j=1:maxiter
        xQEIF_previous=xQEIF;
        xQEIF=ProxQgammaiota(k,gamma,1/t,xQEIF-t*(AtA*xQEIF-A'*b));
        if norm(xQEIF-xQEIF_previous)<levelstoppingcriteria
            break
        else
            convergedQEIF=convergedQEIF+1;
        end
    end
    if convergedQEIF == maxiter,
         disp('Not enough iterations QEl0')
     end

%     figure(3);
%     plot(abs(xQEIF));title('Quadratic envelope of indicator functional');axis([1 n 0 1.2*norm(xo,inf)]);hold on; plot(abs(xo),'r');legend({'oracle solution','reconstruction for different \mu'});hold off
%     pause(0.3)
id=find(abs(xQEIF)>10^(-7));supp(id)=1;
if norm(truesupp-supp,1)<smf_QEIF,
    smf_QEIF=norm(truesupp-supp,1);
    lambda_smf_QEIF=gamma;
end
if norm(xgt-xQEIF)<dgt_QEIF,
    dgt_QEIF=norm(xgt-xQEIF);
    lambda_dgt_QEIF=gamma;
end
 if norm(xo-xQEIF)<100*levelstoppingcriteria,
     count=count+1;
     gamma_foundoracle_QEIF(count)=gamma;
 end


%we dont have a stopping criterion here since anyways we test over very few
%values of gamma
end

smf_QEIF_numtrial_nl(l,iter_nl)=smf_QEIF;
dgt_QEIF_numtrial_nl(l,iter_nl)=dgt_QEIF;
gamma_smf_QEIF_numtrial_nl(l,iter_nl)=lambda_smf_QEIF;
gamma_dgt_QEIF_numtrial_nl(l,iter_nl)=lambda_dgt_QEIF;
gamma_foundoracle_QEIFtot=[gamma_foundoracle_QEIFtot;gamma_foundoracle_QEIF];

end
toc
end

