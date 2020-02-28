%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MBiTe LAB, University of Nebraska Lincoln
% Second run this code (2)
% Convert MS_Matrix.mat,FBAorder and ReducedReactionLabel into single
% MS_2states.txt file where you can load in excel to work as
% MS_2states.xlsx
% Matlab code Developer: Zahmeeth Sakkaff
% Date 06/08/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
load('MS_Matrix.mat');
load('FBAorder');
load('ReducedReactionLabel');

twostates = [FBAOrd128list;num2cell(MS_Matrix)];
s1 = size(twostates);
Original_MS_2states = cell(s1(1),s1(2)+1);
Original_MS_2states(1,2:s1(2)+1) = twostates(1,:);

for i=2:s1(1)
Original_MS_2states(i,:) = [ReducedReactionLabel(i-1),twostates(i,:)];
end

fid = fopen('MS_2states.txt','w');
s2 = size(Original_MS_2states);
for columns = 1:size(Original_MS_2states,2)
      fprintf(fid,'%s\t',Original_MS_2states{1,columns});
end
fprintf(fid,'\n');

for rows = 2:s2(1)   
   fprintf(fid,'%s\t',Original_MS_2states{rows,1});
      for columns = 2:s2(2)
        fprintf(fid,'%d\t',Original_MS_2states{rows,columns});
      end
    fprintf(fid,'\n');
end
fclose(fid);
msg1 = 'MS_2states.txt is created to use in excel and run MutualInfo_MS_2states.m file to calculate';
msg2 = 'Mutual information of Stage I (intracellular) with respect to 7 input compounds .... ';
disp(msg1);
disp(msg2);


    
