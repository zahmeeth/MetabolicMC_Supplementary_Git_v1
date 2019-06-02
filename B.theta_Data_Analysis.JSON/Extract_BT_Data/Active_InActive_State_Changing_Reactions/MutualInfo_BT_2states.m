%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MBiTe LAB, University of Nebraska Lincoln
% Third run this code (3)
% Compute the stage I (intracecullar mutual information) for 7 compund inputs 
% Matlab code Developer: Zahmeeth Sakkaff
% Date 06/08/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% A = xlsread('FBAsMatrix_Corrected'); %Old excel matrix
% Matrix = A;

% Temp = load('FBAMatrix.mat'); %New FBAMatrix matrix
% Matrix = Temp.Matrix;

Temp = load('BT_Matrix.mat');
Matrix = Temp.BT_Matrix;

s1 = size(Matrix);

% Reduce duplicate rows to obtain the reduced Matrix without duplicates in the rows of the Matrix               
NoDuplicateRows = unique(Matrix,'rows','stable');
s2 = size(NoDuplicateRows);     % 173x121
NumDuplicateRows = s1(1)-s2(1) % 251-173 = 78 eliminated 78 duplicate rows

% Print the NumDuplicateRows
for i=1:s1(1)
CompareRow(i,:) = ismember(Matrix,Matrix(i,:),'rows');
end

A = sum(CompareRow);
M = max(A);
count = 0;
for i = 1:s1(1)
    if(A(1,i)> 1)
        count = count +1;
        IndexRow(count,1) = i;
    end
end
DuplicatesRows = zeros(length(IndexRow),M);
s3 = size(DuplicatesRows);
% 
row = 1;
for i =1:s3(1)
P = find(CompareRow(:,(IndexRow(i,1))));
s4 = size(P);
column=1;
    for k = 1:s4(1)
         DuplicatesRows(row,column)=P(k,1);
        column = column+1;
    end
    row = row+1;
end
DetailsofDuplicateRows = unique(DuplicatesRows,'rows','stable')


% Count duplicate columns in the reduced Matrix without duplicates in the rows of the Matrix      
NoDuplicateColumns = unique(NoDuplicateRows','rows','stable');
s5 = size(NoDuplicateColumns) ;
% 
NumDuplicateColumns = s2(2)-s5(1) % 121-117 = 4
% 
% % Print the NumDuplicateColumns
ColumnData = NoDuplicateRows';
for i=1:s2(2)
CompareColumns(i,:) = ismember(ColumnData,ColumnData(i,:),'rows');
end

A = sum(CompareColumns);
M = max(A);
count = 0;
for i = 1:s2(2)
    if(A(1,i)> 1)
        count = count +1;
        IndexColumn(count,1) = i;
    end
end
DuplicateColumns = zeros(length(IndexColumn),M);
s6 = size(DuplicateColumns);
% 
row = 1;
for i =1:s6(1)
P = find(CompareColumns(:,(IndexColumn(i,1))));
s7 = size(P);
column=1;
    for k = 1:s7(1)
         DuplicateColumns(row,column)=P(k,1);
        column = column+1;
    end
    row = row+1;
end
DetailsofDuplicateColumns = unique(DuplicateColumns,'rows','stable')
% % Ans1 = log2(s2(2)) %6.9189
% % Ans2 =((count/2) * ((count)/s2(2)) * (2*(2*0.5))) %0.5289

%  UpperBoundMutualInfor = log2(s2(2)) -((count/2) * ((count)/s2(2)) * (2*(2*0.5))) % New = 6.3899, old = %5.7288
% % % MutualInfor = log2(121) -(6 * (12/121) * (2*(2*0.5)))  %5.7288 (Old)
% % % MutualInfor = log2(121) -(4 * (8/121) * (2*(2*0.5)))  %6.3899  (New)
% MutualInfor = log2(128) -(64 * ((128)/128) * (2*(2*0.5)))  %6.3899  (New)
s8 = size(DetailsofDuplicateColumns);
for i = 1:s8(1)
NoE(i,1) = nnz(DetailsofDuplicateColumns(i,:));%NoE - number of elements in each row in DetailsofDuplicateColumns matrix
end

[CountEle,Ele]=hist(NoE,unique(NoE'))

s9 = size(Ele);
temp = 1;
for i = 1:s9(2)
    Val=0;
    for j = 1:Ele(1,i)
        Val = Val + ((1/Ele(1,i))*log2(1/Ele(1,i)));
    end
    H_XgivenY(temp,1) = CountEle(1,i) * ((Ele(1,i)/s2(2)) * Val);
    temp = temp +1;
end

msg1 = 'Entropy value as in page 12 in the manusript...';  
disp(msg1)
ToTH_XgivenY = -1 *(sum(H_XgivenY))
msg2 = 'Mutual Information value as in page 12 equation 24 in the manusript...';  
disp(msg2)
UpperBoundMutualInfor = log2(s2(2))- ToTH_XgivenY


clearvars -except Data Matrix  DetailsofDuplicateColumns NoE