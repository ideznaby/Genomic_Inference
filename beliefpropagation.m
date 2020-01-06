function [ accuracy, entropyavg, error, foundsnps ] = beliefpropagation( data, dataknown, Fr,FrChild, MAF, familynum, genenum, numofloops, markovorder,hiddenmember, trans, haplotype, genDistance, withphenotypes, phenotypeindices, phenotypevalues, withHMM, knownfamilymembers, familytree )
%BIELEFEPROPAGATION is the function which creates the graph and calculates
%the value of the missing parameters in data which are indicated by -1 and
%then calculates accuracy and average of entropy and outputs them.
%   INPUT data: is the full version of SNP data of 11 members of CEPH/Utah
%   pedigree 1463 family. we use data to calculate the accuracy of the
%   algorithm
%   INPUT dataknown: is the input data in our case dataknown is the SNP data of 11
%   family members of CEPH/Utah pedigree 1463 family and the missing data is indicated as -1 in it.
%   our goal here is to find values indicated by -1.
%   INPUT Fr: is the probabilities of Mendelian's inheritence law
%   INPUT FrChild: is the probabilities of Mendelian's inheritence law when the
%   child is given
%   INPUT MAF: major allele frequency, shows at each SNP what is the
%   probability of each allele happening
%   INPUT familynum: is the number of familly members which in our case is 11
%   genenum is the number of SNPs for any individual.
%   INPUT numofloops: indicate how many loops shoud loopy belief propagation go
%   markovorder shows which markovorder should the algorithm calculate for
%   g nodes.
%   INPUT hiddenmember: if hiddenmember is not 0 this function will only
%   calculates the accuracy and entropy corresponding to the indicated
%   familly member by this variable
%   INPUT trans: Calculated transition probabilities should passed to
%   method with this argument the probabilities are calculated using
%   probe-calculate method.
%   INPUT haplotype: The haplotype data of the population
%   INPUT genDistance: The gene distance data 
%   INPUT withphenotypes: a boolean value indicating whether or not to use phenotype data
%   INPUT phenotypeindices: The SNP indices of phenotypes we know about family members
%   INPUT phenotypevalues: Values of the SNPs which we know their indices (this should be a 3 dimentional matrix)
%   where the first index is for 3 values of the SNP (0,1,2) the second one is for each phenotype
%   so the number here should be equal to the length of phenotypeindices and the last one is for each family member.
%   INPUT withHMM: boolean indicating if the program should use recombination model or Markov chain (if it is true the model will use recombination model) 
%   INPUT knownfamilymembers: Indicates the family members who should be
%   included in the graph so the target should always be here.
x = {};
xe = {};
f = {};
g = {};
id = 1;
%this indicates that our belief propagation has loop inside.
global loopy;
loopy = 1;
%creates familialnodes.
f= cell(size(familytree,1),genenum);
for i=1:size(familytree,1)
    for j=1:genenum
        if (any(knownfamilymembers == familytree(i,3)) && any(knownfamilymembers == familytree(i,1))) || (any(knownfamilymembers == familytree(i,1)) && any(knownfamilymembers == familytree(i,2)))% Only if one of the relations are in dataknown create the familial node
            f{i,j} = familialnode(id);
            id = id+1;
            f{i,j}.Fr = Fr;
            f{i,j}.FrChild = FrChild;
            f{i,j}.MAF = MAF;
            f{i,j}.num_of_individuals = familynum;
            f{i,j}.numofSNPs = genenum;
            f{i,j}.knownfamilymembers = knownfamilymembers;
        end
    end
end
xeids = 10000000;
%create g and x nodes and connect them to eachother and f nodes.
for i=1:familynum
    if any(knownfamilymembers == i)
        g{i} = gnode(id);
        g{i}.isHMM = withHMM;
        g{i}.trans = trans;
        g{i}.haplotype = haplotype;
        g{i}.genDistance = genDistance;
        g{i}.markovorder = markovorder;
        g{i}.numofSNPs = genenum;
        g{i}.knownfamilymembers = knownfamilymembers;
        id =id+1;
    end
    for j=1:genenum
        x{i,j} = evident_node(id);
        id = id+1;
        xe{i,j} = equ_node(xeids);
        xe{i,j}.MAF = MAF(j,:);
        xe{i,j}.knownfamilymembers = knownfamilymembers;
        xe{i,j}.withHMM = withHMM;
        if any(knownfamilymembers == i)
            connect(xe{i,j},g{i});
        end
        connect(x{i,j},xe{i,j});
        %connecting every node to appropriate familial node
        for k=1:size(familytree,1)
            if ~isempty(f{k,j})
                if i == familytree(k,1) || i == familytree(k,2) || i == familytree(k,3)
                    connect(xe{i,j},f{k,j});
                    if i == familytree(k,3)
                        f{k,j}.childnodeid = xeids;
                    end
                end
            end
        end
        xeids = xeids + 1;
    end
end
%adding phenotypes:
if withphenotypes
    ph = cell(size(dataknown,1),length(phenotypeindices));% we will assume we know the same phenotypes for all the family members
    for i = knownfamilymembers
        for j=1:length(phenotypeindices)
            id = id + 1;
            ph{i,j} = evident_node(id);
            connect(ph{i,j},xe{i,phenotypeindices(j)});
        end
    end
end
unknown = [];
%setting the initial messages of the nodes.
for i=1:familynum
    for j=1:size(dataknown,2)
        if dataknown(i,j) == 0
            x{i,j}.setup_init_msg([1 0 0]);
        end
        if dataknown(i,j) == 1
            x{i,j}.setup_init_msg([0 1 0]);
        end
        if dataknown(i,j) == 2
            x{i,j}.setup_init_msg([0 0 1]);
        end
        if dataknown(i,j) == -1
            %ms = MAF(j,:) + 0.01;
            %x{i,j}.setup_init_msg(ms / sum(ms));
            x{i,j}.setup_init_msg([1 1 1]);
            if(hiddenmember == 0)
                unknown = [unknown;i,j];
            else
                if i == hiddenmember
                    unknown = [unknown;i,j];
                end
            end
        end
    end
end
if withphenotypes
    %Setting the messages of phenotype nodes
    for i=knownfamilymembers
        for j=1:length(phenotypeindices)
            ph{i,j}.setup_init_msg(phenotypevalues(:,j,i)');
        end
    end
end
%update every node every time for number of loops given
% if we increase the numofloops we have more chance the values will
% converge and our accuracy will increase
preventavg([1,2,3]) = 2;
preventropy = zeros(1,size(unknown,1));
loopcout = 0;
prevent = 2;
for k = 1:numofloops
    disp(['loop number : ' , num2str(k)])
    for i = 1:familynum
        for j = 1:size(xe,2)
            xe{i,j}.update_node([1/3 1/3 1/3]);
        end
    end
    for i= 1:familynum
        if any(knownfamilymembers == i)
            g{i}.update_node([1/3 1/3 1/3]);
        end
    end
    for i = 1:size(familytree,1)
        if ~isempty(f{i,1})
            for j = 1:size(f,2)
                f{i,j}.update_node([1/3 1/3 1/3]);
            end
        end
    end
    entropysum = 0;
    for i = 1: size(unknown,1)
        foundsnp = (marginal(x{unknown(i,1),unknown(i,2)},xe{unknown(i,1),unknown(i,2)}));
        sument = 0;
        for j=1:3
            ent = log2(foundsnp(j)) * foundsnp(j);
            if foundsnp(j) == 0
                ent = 0;
                continue;
            end
            sument = sument + ent;
        end
        entropy(i) = sument / log2(3);
        entropysum = entropysum - entropy(i);
    end
    entropyavg = entropysum / size(unknown,1);
    % if entropy does not change we stop the loop
    if abs(prevent - entropyavg) < 0.00000001
        break;
    end
    % if after 10 loops the entropy starts to oscillate for 30 loops we
    % break
    if abs(preventavg(mod(k-2,3)+1) - entropyavg) < 0.00000001 && k > 10
        loopcout = loopcout + 1;
    else
        loopcout = 0;
    end
    if loopcout > 30
        break;
    end
    preventavg(mod(k,3)+1) = entropyavg;
    prevent = entropyavg;
    preventropy = entropy;
end
% find marginals values of each link
disp('************************')
true = 0;
foundsnps = zeros(1,size(unknown,1));
%calculate the accuracy
for i = 1: size(unknown,1)
    [~,foundsnps(i)] = max(marginal(x{unknown(i,1),unknown(i,2)},xe{unknown(i,1),unknown(i,2)}));
    if foundsnps(i)-1 == data(unknown(i,1),unknown(i,2))
        true = true + 1;
    end
end
accuracy = true / size(unknown,1);
%calculate the correctness
allerror = 0;
classes = 0:2;
entropysum = 0;
for i = 1: size(unknown,1)
    foundprob = (marginal(x{unknown(i,1),unknown(i,2)},xe{unknown(i,1),unknown(i,2)}));
    sument = 0;
    sumentMAF = 0;
    for j=1:3
        ent = log2(foundprob(j)) * foundprob(j);
        if foundprob(j) == 0
            ent = 0;
            continue;
        end
        sument = sument + ent;
    end
    entropy(i) = sument / log2(3);
    [~,foundsnp] = max(foundprob);
    errorsum = 0;
    for j = 0:2
        errorsum = errorsum + foundprob(j+1) .* double(abs(int32(j)-data(unknown(i,1),unknown(i,2))));
    end
    allerror = allerror + errorsum;
    entropysum = entropysum - entropy(i);
end
error = allerror / size(unknown,1);
%calculate entropy    
entropyavg = entropysum / size(unknown,1);
end

