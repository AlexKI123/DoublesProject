% assign single-double-triple mutational class to input data on TP53 drivers
% 
% gencode.AA is from the file gencode.txt
%
% input data in 'data_P53.txt' is from Supp Table 3 of Giacomelli, et al
% ('https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0204-y/MediaObjects/41588_2018_204_MOESM6_ESM.xlsx')

for j=1:length(data_P53.codon)
    % extract codon from data
    tri=char(data_P53.codon(j));
    % check single mutation neighbours
    pos1 = gencode.AA(find(gencode.FirstPosition~=tri(1)&gencode.SecondPosition==tri(2)&gencode.ThirdPosition==tri(3)));
    pos2 = gencode.AA(find(gencode.FirstPosition==tri(1)&gencode.SecondPosition~=tri(2)&gencode.ThirdPosition==tri(3)));
    pos3 = gencode.AA(find(gencode.FirstPosition==tri(1)&gencode.SecondPosition==tri(2)&gencode.ThirdPosition~=tri(3)));
    % check tandem mutation neighbours
    pos12 = gencode.AA(find(gencode.FirstPosition==tri(1)&gencode.SecondPosition~=tri(2)&gencode.ThirdPosition~=tri(3)));
    pos23 = gencode.AA(find(gencode.FirstPosition~=tri(1)&gencode.SecondPosition~=tri(2)&gencode.ThirdPosition==tri(3)));
    if (sum(find(pos1==data_P53.AA_variant(j)))>0)||(sum(find(pos2==data_P53.AA_variant(j)))>0)||(sum(find(pos3==data_P53.AA_variant(j)))>0)
        data_P53.mutation(j)="single";
    elseif (sum(find(pos12==data_P53.AA_variant(j)))>0)||(sum(find(pos23==data_P53.AA_variant(j)))>0)
        data_P53.mutation(j)="tandem";
    else
        data_P53.mutation(j)="triplet";
    end
end
