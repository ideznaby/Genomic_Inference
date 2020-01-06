clear
clc
%% adding the path of the codes
addpath('SumProductLab_R3.00');%The codes for the graph nodes is in this folder it is taken from SumProductLab from url: https://www.mathworks.com/matlabcentral/fileexchange/26607-sumproductlab-for-factor-graphs
addpath('HiddenMarkov'); %This folder contains the code for Hidden Markov Model (recombination model)
addpath('Observable Markov');% This folder contains the code for Markov Chain part
addpath('Data');%This folder contains the data
%% Loading data
load('UTAH_Pedigree_1463_data.mat'); % Load the family data
% in this data the rows are in the following order:
% 1- The ID of SNPs 2- Data of P5 3- P6 4- C7 5- C8 6- C9 7- C10 8- C11 9- GP1 10- GP2 11- GP3 12- GP4
familySNPs = Utah_Family_Data;
IDs = familySNPs(1,:);
familytree = ReadFamilyTree('UTAH_Pedigree_1463_Family_Tree.txt');
%% load the calculated major allele frequency over all of the data.
%[MAF ib1] = ReadMAF('extractedMAF.txt',IDs, true);
load('Genome_1000_maf_intersected.mat');
MAF = Genome_1000_maf_intersected;
[IDs,ia,ib] = intersect(MAF(:,1),IDs);
MAF = MAF(ia,:);
familySNPs = familySNPs(:,ib);
% loading the data of population
[populationdata indices IDs] = ReadPopulation('1000genomepopdata_intersected.txt',IDs);
familySNPs = familySNPs(:,indices);
%% intersect MAF with the family data so we only have the data which we have MAF of them
[IDs ia ib] = intersect(MAF(:,1),IDs);
MAF = MAF(ia,:);
familySNPs = familySNPs(:,ib);
%% Defining Some parameters for the model such as Mendel's inheritance probabilites, number of family members, number of SNPs, number of loops and markov order
% Mendel's Inheritance probabilities ( Figure 1 in Paper)
Fr(1,1,:) = [1 0 0];
Fr(1,2,:) = [0.5 0.5 0];
Fr(1,3,:) = [0 1 0];
Fr(2,1,:) = [0.5 0.5 0];
Fr(2,2,:) = [0.25 0.5 0.25];
Fr(2,3,:) = [0 0.5 0.5];
Fr(3,1,:) = [0 1 0];
Fr(3,2,:) = [0 0.5 0.5];
Fr(3,3,:) = [0 0 1];
% FR for the child
FrChild(1,1,:) = [0.5 0.5 0];
FrChild(1,2,:) = [0.5 0.5 0];
FrChild(1,3,:) = [0 0 0];
FrChild(2,1,:) = [0 0.5 0.5];
FrChild(2,2,:) = [0.33 0.33 0.33];
FrChild(2,3,:) = [0.5 0.5 0];
FrChild(3,1,:) = [0 0 0];
FrChild(3,2,:) = [0 0.5 0.5];
FrChild(3,3,:) = [0 0.5 0.5];
% number of family members we will use
familynum = 11;
% number of SNPs we will use
genenum = 100;
% number of loops for our belief propagation
numofloops = 50;
% order of markov chain
markovorder = 3;
%% Pick the first 100 SNPs
MAF = MAF(1:genenum,2:end);
intersectedFamilySNPs = familySNPs(:,1:genenum);
populationdata = populationdata(1:genenum,:);
data = intersectedFamilySNPs(2:end,:);
%% There is no need for loading genDistances and haplotypes for now because in first example we are not using recombination model.
%genDistance = genDistance(1:genenum,:);
genDistance = [];

%haplotype = haplotype(1:genenum,:);
haplotype = [];
%% EXAMPLE1: Setting parameters for a sample run of the algorithm with markov chain with order 3
withHMM = false; % should the program use recombination model?
usephenotypes = true; % should the program use phenotypes?
markovorder = 3; % the markov chain order
trans = probe_calculate(markovorder,populationdata); % calculate the transitional probabilities for markov chain using the population data
trans = trans + 0.01; % adding some noise to these transitions ( they are not probabilities anymore but the algorithm manages them so we should not worry about it now)
dataknown = data; % dataknown is our input data for places that we want to predict we put -1 and the rest will be the real values
dataknown(1,:) = removeatrandom(data(1,:),0.5); % put -1 for random 50 percent of the SNPs the first member(P5 in this case) 
knownfamilymembers = [1,2,3,4,5,6,7,8,9,10,11]; % indicate which family members should be included in the graph
hiddenmember = 1; % the member we want to predict its SNPs in this case P5 note that the program predicts all -1s in the dataknown however, calculates outputs only on this member ( put 0 if you want to calculate over all members)

% Defining Phenotypes
phenotypeindices = [65 52 34];% indices for some test phenotypes
phenotypevalues = zeros(3,3,11);%the first index is for 3 values of the SNP (0,1,2) the second one is for each phenotype so the number here should be equal to the length of phenotypeindices and the last one is for each family member
for i=1:11 %Simulating the phenotypes from the data we know about each family member
    for j=1:length(phenotypeindices)
        phenovalue = ([0 1 2] == data(i,phenotypeindices(j)));
        phenovalue = (phenovalue + 0.1) / 1.3;%add noise to phenotype values but make sure they are still probability
        phenotypevalues(:,j,i) = phenovalue;
    end
end
% Running the algorithm with given parameters
[Accuracy,Entropy,Error,~] = beliefpropagation( data, dataknown, Fr,FrChild, MAF,familynum, genenum, numofloops, markovorder,hiddenmember, trans , haplotype, genDistance, usephenotypes, phenotypeindices, phenotypevalues, withHMM, knownfamilymembers, familytree)




%% EXAMPLE2: Reading files and running the algorithm with recombination model
% Reading the data again because the previous data is intersected with
% population data
load('UTAH_Pedigree_1463_data.mat'); % Load the family data
% in this data the rows are in the following order:
% 1- The ID of SNPs 2- Data of P5 3- P6 4- C7 5- C8 6- C9 7- C10 8- C11 9- GP1 10- GP2 11- GP3 12- GP4
familySNPs = Utah_Family_Data;
IDs = familySNPs(1,:);
%%
% reading genDistance and haplotypes
fid = fopen('small_genetic_map_chr22_combined_b36.csv'); % Load the genetic distance data for Recombination Model
genDistance = cell2mat(textscan(fid,'%f\t%f'));
file = fopen('small_CEU.chr22.csv'); %Read Haplotype data for recombination model
format = repmat('%d ',1,235);
haplotype = cell2mat(textscan(file,format)); 
%% Intersect the familySNPs data IDs with genDistance and haplotypes so we only left with the SNPs we have all the data for
[IDs,ia,ib] = intersect(IDs,genDistance(:,1));
genDistance = genDistance(ib,:);
familySNPs = familySNPs(:,ia);
[IDs,ia, ib] = intersect(IDs,haplotype(:,1));
haplotype = haplotype(ib,:);
familySNPs = familySNPs(:,ia);
%% Loading and intersecting the new MAF files which is intersected with genDistance and Haplotypes
[MAF ib] = ReadMAF('extractedMAF.txt',IDs, true);
familySNPs = familySNPs(:,ib);
%% Picking the first 100 haplotypes and genDistances
genDistance = genDistance(1:genenum,:);
haplotype = haplotype(1:genenum,:);
intersectedFamilySNPs = familySNPs(:,1:genenum);
data = intersectedFamilySNPs(2:end,:);

%% Running the algorithm with Recombination Model
withHMM = true; % should the program use recombination model?
usephenotypes = false; % should the program use phenotypes?
markovorder = []; % the markov chain order
populationdata = []; % We do not need population data for recombination model
trans = [];
dataknown = data; % dataknown is our input data for places that we want to predict we put -1 and the rest will be the real values
dataknown(1,:) = removeatrandom(data(1,:),0.5); % put -1 for random 50 percent of the SNPs the first member(P5 in this case) 
knownfamilymembers = [1,2,3,4,5,6,7,8,9,10,11]; % indicate which family members should be included in the graph
hiddenmember = 1; % the member we want to predict its SNPs in this case P5 note that the program predicts all -1s in the dataknown however, calculates outputs only on this member ( put 0 if you want to calculate over all members)

% Defining Phenotypes
phenotypeindices = [];% Since usephenotypes is false we do not need these
phenotypevalues = [];
% Running the algorithm with given parameters
[Accuracy,Entropy,Error,~] = beliefpropagation( data, dataknown, Fr,FrChild, MAF,familynum, genenum, numofloops, markovorder,hiddenmember, trans , haplotype, genDistance, usephenotypes, phenotypeindices, phenotypevalues, withHMM, knownfamilymembers, familytree)