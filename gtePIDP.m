function maxpop = gtePIDP(M,E,D,I,k,lambda,p_threshold,exclusion)

%
% This function uses genetic algorithm to solve the integrative model.
%
% M : (sample) x (mutation_gene)         mutation matrix
% E : (expression_gene) x (sample)       expression matrix
% D : (expression_gene) x (TF,miRNA)     TF and miRNA targeting matrix
% I : (TF,miRNA) x (mutation_gene)       the impact matrix of somatic alterations on TFs and miRNAs
%
% k : number of genes to be identified     needs to be specified by the user
%
% lambda : the parameter in the new measure
% p_threshold : i.e., the parameter \beta in the model
% exclusion : the genes which are to be excluded from the results
%
% maxpop : output matrix with k+2 columns.
%          every row indicates the information of one solution; 
%          for row i, maxpop(i,1:k) records the selected genes;
%          maxpop(i,k+1) records the weight of this gene set;
%          maxpop(i,k+2) records the significance level. 
 
[m,n]=size(M);
[m_e,n_e]=size(E);   % Actually, n_e=m
[m_t,n_t]=size(D);   % Actually, m_t=n_e
 
if n<100
    popsize=n;   % the number of individual in one population
elseif n>=100
    popsize=100*(floor(n/100));
end
nger=1000;                  % maximal number of iterations, can be adjusted by users
 
if ~isempty(exclusion)
    M(:,exclusion)=0;
end
 
object_value=zeros(1,nger);
pop=zeros(popsize*2,k+6);   % population matrix, every row record an individual
 
% Generate initial population
for i=1:popsize*2
    temp=randperm(n);
    pop(i,1:k)=temp(1:k);
    pop(i,k+1)=fit_M_E(M,E,D,I,pop(i,1:k),lambda,p_threshold);
end
[~,J]=sort(pop(:,k+1),'descend');
pop=pop(J,:);              % sort the population in descending order
 
R=0; % record whether current optimal solution is trapped in one local optimal solution
i=1;
while(R<10 && i<=nger) % termination
    fit_vector=pop(:,k+1);
    temp_maxweight=max(fit_vector);
    j=1;
    while(j<=popsize) % generate one new individual
        [index1 index2]=select_order_fitness(fit_vector);
        pop(popsize+j,1:k)=crossover(pop(index1,1:k),pop(index2,1:k),n); % crossover
        pop(popsize+j,1:k)=mutation_SA_M_E(M,E,D,I,pop(popsize+j,1:k),n,1,lambda,p_threshold);  % mutation
        pop(popsize+j,k+1)=fit_M_E(M,E,D,I,pop(popsize+j,1:k),lambda,p_threshold);
        j=j+1;
     end
    [~,J]=sort(pop(:,k+1),'descend');
    pop=pop(J,:); % sort the population in descending order, the maximal n individual are transfered to next generation
    object_value(i)=pop(1,k+1);
    
    % use local search to improve the current optimal solution when it is trapped in one local optimal solution
    if R==2
        for j=1:popsize*0.05
            pop(j,1:k)=mutation_SA_M_E(M,E,D,I,pop(j,1:k),n,sqrt(n),lambda,p_threshold);
            pop(j,k+1)=fit_M_E(M,E,D,I,pop(j,1:k),lambda,p_threshold);
        end
    end 
    
    if R==5
        for j=1:popsize*0.01
            pop(j,1:k)=mutation_SA_M_E(M,E,D,I,pop(j,1:k),n,n,lambda,p_threshold);
            pop(j,k+1)=fit_M_E(M,E,D,I,pop(j,1:k),lambda,p_threshold);
        end
    end
    
    maxweight=max(pop(:,k+1));
 
    if maxweight==temp_maxweight
        R=R+1;
    else
        R=0;
    end
    
    i=i+1;
end
[~,J]=sort(pop(:,k+1),'descend');
pop=pop(J,:);
maxpop=pop(1:popsize,:);
 
% delete the repetitive solution
[m,~]=size(maxpop);
i=1;
while(i<m)
    index=zeros(1,m);
    for j=i+1:m
        if all(maxpop(j,:)==maxpop(i,:))
            index(j)=1;
        end
    end
    index=logical(index);
    maxpop(index,:)=[];
    [m,~]=size(maxpop);
    i=i+1;
end
 
% significance test
[m,~]=size(maxpop);
for i=1:m
    maxpop(i,k+2)=fit_M(M,maxpop(i,1:k));
    [Result1,Result2]=impacted_acctivity(M,E,D,I,maxpop(i,1:k),p_threshold);
    maxpop(i,k+3)=Result1+Result2;
    maxpop(i,k+4)=Result1;
    maxpop(i,k+5)=Result2;
    [num_e,num_u,vec_e,vec_u]=impacted_acctivity_num(M,E,D,I,maxpop(i,1:k),p_threshold);
    maxpop(i,k+6)=num_e;
    maxpop(i,k+7)=num_u;
    maxpop(i,k+8)=significance_M_E(M,E,D,I,maxpop(i,1:k),lambda,p_threshold);
    leng_e=length(vec_e);
    leng_u=length(vec_u);
    maxpop(i,k+9:k+8+leng_e)=vec_e;
    maxpop(i,k+9+leng_e:k+8+leng_e+leng_u)=vec_u;
end
 
 
%% The fuctions used in the above main function
function [Result1,Result2]=impacted_acctivity(M,E,D,I,x,p_threshold)  % calculate the impacted_acctivity
if ~all(x)
    x=find(x==1);
end
k=length(x);
 
M1=M;
I1=I;
 
M1(:,x)=[];
I1(:,x)=[];
 
E1=D*I1*M1';
P_e=mattest(E1,E);
 
A=M*I';
A1=M1*I1';
P_u=mattest(A1',A');
 
r_e=length(P_e);
S_e=[];
[val_e,loc_e]=sort(P_e);
for i=1:r_e
    if val_e(i) < p_threshold
        S_e=[S_e val_e(i)];
    else
        break
    end
end
 
r_u=length(P_u);
S_u=[];
[val_u,loc_u]=sort(P_u);
for i=1:r_u
    if val_u(i) < p_threshold
        S_u=[S_u val_u(i)];
    else
        break
    end
end
 
l_e=length(S_e);
l_u=length(S_u);
 
z_e=0;
for i=1:l_e
    z_e=z_e+norminv(1-S_e(i),0,1);
end
Result1=z_e/sqrt(l_e);
 
z_u=0;
for i=1:l_u
    z_u=z_u+norminv(1-S_u(i),0,1);
end
Result2=z_u/sqrt(l_u);
 
 
function [num_e,num_u,vec_e,vec_u]=impacted_acctivity_num(M,E,D,I,x,p_threshold)  % calculate the numbers of the impacted_acctivity
 
if ~all(x)
    x=find(x==1);
end
k=length(x);
 
M1=M;
I1=I;
 
M1(:,x)=[];
I1(:,x)=[];
 
E1=D*I1*M1';
P_e=mattest(E1,E);
 
A=M*I';
A1=M1*I1';
P_u=mattest(A1',A');
 
r_e=length(P_e);
S_e=[];
vec_e=[];
[val_e,loc_e]=sort(P_e);
for i=1:r_e
    if val_e(i) < p_threshold
        S_e=[S_e val_e(i)];
        vec_e=[vec_e loc_e(i)];
    else
        break
    end
end
 
r_u=length(P_u);
S_u=[];
vec_u=[];
[val_u,loc_u]=sort(P_u);
for i=1:r_u
    if val_u(i) < p_threshold
        S_u=[S_u val_u(i)];
        vec_u=[vec_u loc_u(i)];
    else
        break
    end
end
 
num_e=length(S_e);
num_u=length(S_u);
 
 
function f=fit_M_E(M,E,D,I,x,lambda,p_threshold)   % fitness calculation function
 
k=length(x);
[~,n]=size(M);
temp=zeros(1,n);
temp(x)=1;
x=temp;
index = x==1;
M_index=M(:,index);
M_index_sum=sum(M_index,2);
[Result1,Result2]=impacted_acctivity(M,E,D,I,x,p_threshold);
f=(2*sum(M_index_sum>0)-sum(M_index_sum))+lambda*(Result1+Result2);
 
 
function f=fit_M(M,x)   % mutation fitness calculation function
 
k=length(x);
[~,n]=size(M);
temp=zeros(1,n);
temp(x)=1;
x=temp;
index = x==1;
M_index=M(:,index);
M_index_sum=sum(M_index,2);
f=2*sum(M_index_sum>0)-sum(M_index_sum);
 
 
function [index1 index2]=select_order_fitness(fit_vector) % select two individuals based on order fitness
n=length(fit_vector);
[~,J]=sort(fit_vector);
p=zeros(1,n);
for i=1:n
    p(J(i))=2*i/(n*(n+1));
end
pp=cumsum(p); 
random_data=rand(1,1);
temp=find(pp>=random_data);
index1=temp(1);
random_data=rand(1,1);
temp=find(pp>=random_data);
index2=temp(1);
 
 
function newpop=crossover(parent1,parent2,n)  % crossover function
temp=zeros(1,n);
temp(parent1)=1;
parent1=temp;
temp=zeros(1,n);
temp(parent2)=1;
parent2=temp;
k=sum(parent1);
newpop=zeros(1,n);
index=(parent1+parent2==2);
newpop(index)=1;
parent1(index)=0;
parent2(index)=0;
temp=find(parent1+parent2==1);
index=randperm(sum(parent1+parent2==1));
newpop(temp(index(1:(k-sum(newpop)))))=1;
newpop=find(newpop==1);
 
 
function m_x=mutation(x,n) % mutation function
temp=zeros(1,n);
temp(x)=1;
x=temp;
k=sum(x);
index1=randi(n,1);
while(x(index1)==1)
    index1=randi(n,1);
end
x_nonzero=find(x==1);
index2=x_nonzero(randi(k,1));
m_x=x;
m_x(index1)=1;
m_x(index2)=0;
m_x=find(m_x==1);
 
 
function m_x=mutation_SA_M_E(M,E,D,I,pop_i,n,N,lambda,p_threshold)   % local search
for i=1:N
    pop_j=mutation(pop_i,n);
    if fit_M_E(M,E,D,I,pop_j,lambda,p_threshold)>=fit_M_E(M,E,D,I,pop_i,lambda,p_threshold)
        pop_i=pop_j;
    end
end
m_x=pop_i;
 
 
function p=significance_M_E(M,E,D,I,subset,lambda,p_threshold)       % significance test function
[m,~]=size(M);
w=zeros(1,10);
n=length(subset);
for j=1:10
    M_temp=M;
    E_temp=E;
    M_temp(:,subset)=0;
    for i=1:n
        temp=sum(M(:,subset(i)));
        index = randperm(m,temp);
        M_temp(index,subset(i))=1;
    end
    w(j)=fit_M_E(M_temp,E,D,I,subset,lambda,p_threshold);
end
p=sum(w>=fit_M_E(M,E,D,I,subset,lambda,p_threshold))/100;
