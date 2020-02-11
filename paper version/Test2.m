%This script was used to generate the data in Test 2 of the paper
%Perfect support recovery with LASSO. The script generats simulations over
%3 different parameters, m (amount of rows), r (relative amount of non-zeros in
%ground truth) and numtrial (which simply repeats the experiment for same
%value of m and r).


delta_lower_m_r_numtrial=[];%This parameter will store values of 1-sigma_k(A'A) for k=100*r/m 
delta_upper_m_r_numtrial=[];%This parameter will store values of sigma_1(A'A)-1 for k=100*r/m 
proj_m_r_numtrial=[];%Stores values of the norm of the projection of a random normalized vector onto the range of A 
%(which is of relevance for formula (25) in the paper)



sizes=[(1:10) (20:10:100) (200:100:400)];

for ind=1:22
    m=50*sizes(ind)
    delta_lower_r_numtrial=[];
    delta_upper_r_numtrial=[];
    proj_r_numtrial=[];
for k=[1:10]*sizes(ind);%since m is 50 times sizes, this means r=2,4... etc until 20%    
    delta_lower_numtrial=[];
    delta_upper_numtrial=[];
    proj_numtrial=[];
for numtrial=1:20
    A=randn(m,k)+i*randn(m,k);
    A=A./vecnorm(A);
    e=randn(m,1)+i*randn(m,1);
    e=e/norm(e);
    B=A'*A;%we form this matrix once for faster evaluation
    delta_lower_numtrial=[delta_lower_numtrial 1-(min(svd(B)))];
    delta_upper_numtrial=[delta_upper_numtrial (max(svd(B)))-1];
    proj_numtrial=[proj_numtrial norm(A*inv(B)*A'*e)];%the formula is well known to yield the projection of e onto range of A
end
    delta_lower_r_numtrial=[delta_lower_r_numtrial; delta_lower_numtrial];
    delta_upper_r_numtrial=[delta_upper_r_numtrial; delta_upper_numtrial];
    proj_r_numtrial=[proj_r_numtrial; proj_numtrial];
end
    delta_lower_m_r_numtrial(ind,:,:)=delta_lower_r_numtrial;
    delta_upper_m_r_numtrial(ind,:,:)=delta_upper_r_numtrial;
    proj_m_r_numtrial(ind,:,:)=proj_r_numtrial;    
end


