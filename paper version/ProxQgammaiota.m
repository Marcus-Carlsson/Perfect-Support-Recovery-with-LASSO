
function x = ProxQgammaiota(k,gamma,rho,y) %k is the cardinality, r=rho and g=gamma, y the point in which we compute prox
%which should be a column vector.

r=abs(y);%sets up polar coordinates
theta=exp(i*angle(y));
n=length(y);

[rsort,id]=sort(r,'descend');
idinv(id)=[1:n];
rnew=[rsort(1:k);rho/gamma*rsort(k+1:n)];%this vector will be the new radii, to be updated before the end

%The following steps follow those described in our paper Carlsson, Marcus, and Daniele Gerosa. 
%"On phase retrieval via matrix completion and the estimation of low rank PSD matrices." Inverse Problems (2019).

%First we compute the proximal operator of the S transform 

%1st step
if (rho/gamma)*rsort(k+1)<rsort(k)
    x=rnew;
    x=x(idinv);x=x.*theta;

else    
    %2nd step: We compute the indices called l* and j* in the paper; we call them just l and j.
    temp=find(rnew<=rnew(k+1));j=min(temp);
    temp=find(rnew>=rnew(k));l=max(temp);
    
    z=sort(rnew(j:l),'descend'); %3rd step. We compute the vector z.
    for m=1:length(z)-1
        s=(z(m)+z(m+1))/2;
        %4th step. We compute the indices called l and j in the paper. We call them l1 and j1.
        temp=find(rnew<=s);j1=min(temp);
        temp=find(rnew>=s);l1=max(temp);
        sI=(rho*sum(rsort(j1:l1)))/((k+1-j1)*rho + (l1-k)*gamma);%5th step
        if (sI>=z(m+1) && z(m)>=sI)%6th step
            x=[max(rnew(1:k),sI);min(rnew(k+1:n),sI)];
            x=x(idinv);x=x.*theta;
            break
        end
    end    
end

%Final step taks prox of S to prox of Q, see formula (20) of the paper
    x=(rho*y-gamma*x)/(rho-gamma);


