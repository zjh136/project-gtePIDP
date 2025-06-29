% We use the package AR_run from the following paper to compute the impact matrix I of somatic alterations on TFs and miRNAs in our EAR model:
% Osmanbeyoglu,H.U., et al. Nat. Commun., 8, 14249 (2017).


% D_1 : (expression_gene) x (TF)         TF targeting matrix
% D_2 : (expression_gene) x (miRNA)      miRNA targeting matrix
% M : (sample) x (mutation_gene)         mutation matrix
% E : (expression_gene) x (sample)       expression matrix

% I : (TF,miRNA) x (mutation_gene)       the impact matrix of somatic alterations on TFs and miRNAs


D=[D_1 D_2];

lambda=0.3;
rsL2=0;
spectrumA=1;
spectrumB=1;
old_version=1;

model = ar_train(D,M,E,lambda, rsL2, spectrumA, spectrumB, old_version);

I = ar_model2w(model);
