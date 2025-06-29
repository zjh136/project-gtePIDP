
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


load M;
load E;
load D;
load I;

lambda=0.3;
p_threshold=0.0001;

exclusion=[];

    
for k=2:10  
    
    maxpop = gtePIDP(A,E,T,L,k,lambda,p_threshold,exclusion);

    eval(['save data',num2str(k)]);
    
end
