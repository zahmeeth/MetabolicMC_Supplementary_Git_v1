%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MBiTe LAB, University of Nebraska Lincoln
% Start with this code (1)
% Construct a MS_Matrix,FBAorder and ReducedReactionLabel from KBase FBA Comparison data file 
% CompareFBA_JmmolMS.json where we can use as input for function FBAComparison.m
% Matlab code Developer: Zahmeeth Sakkaff
% Date 06/08/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% Data = loadjson('CompareFBA_JmmolMS.json');
% save('CompareFBA_JmmolMS.mat','Data');
load('CompareFBA_JmmolMS.mat');
 
s1 = size(Data.reactions);

for i = 1:s1(2)
  Temp.NoFBAs(i) = length(fieldnames(Data.reactions{1,i}.reaction_fluxes));
  Temp.NoStates(i)= length(fieldnames(Data.reactions{1,i}.state_conservation));
  Temp.MCS{i} = char(Data.reactions{1,i}.most_common_state);
end
FBASolu = [num2cell((1:s1(2))'),num2cell((Temp.NoFBAs)'),num2cell((Temp.NoStates)'),cellstr(Temp.MCS)'];

count =1;
%Extract All FBA's less than 100% - All the FBAs take only single state
for i = 1:s1(2)
    if (cell2mat(FBASolu(i,3))>1)
        ReducedFBASolu(count,:)= FBASolu(i,:);
        count = count+1;
    end
end
s2=size(ReducedFBASolu); %556
AllStatus = char(ReducedFBASolu{:,4});
%AllDifStates = 'FOR','IA','NA','REV'
AllDifStates = unique(cellstr(AllStatus),'sorted');

for i = 1:s2(1)
CountStates(i) = ReducedFBASolu{i,3};
end
UniStates = unique(CountStates,'sorted');
% 
MS_Matrix = zeros(s2(1),length(Data.fbas));
countFBAs = 0;
%Total FBAs (127 and less than 127 with 'FOR','IA','NA','REV'status) =556
for i=1:s2(1)
    ReducedIDX(i,1) = ReducedFBASolu{i,1};
    if(ReducedFBASolu{i,2}==length(Data.fbas))%if all the FBAs are ==127 
            countFBAs = countFBAs+1; %with 'FOR','IA','NA','REV'status =433
            P = fieldnames(Data.reactions{1,ReducedFBASolu{i,1}}.reaction_fluxes);
            for j=1:length(P)
            S = strcat('Data.reactions{1,',num2str(ReducedFBASolu{i,1}),'}.reaction_fluxes','.',P{j,1}) ;
            Q = eval(S);
             if(strcmp(Q(1),'FOR'))
                 R(1,j)=1;
                    elseif(strcmp(Q(1),'IA'))
                      R(1,j)=0;
                         elseif(strcmp(Q(1),'NA'))
                           R(1,j)=0;
                                elseif(strcmp(Q(1),'REV'))
                                 R(1,j)=1;
             end                      
            end
       MS_Matrix(i,:) = R;
    end
     
end

%print the FBA labels
Convert = char(P);
for i= 1:length(Data.fbas)
FBAOrd127(i,:) = Convert(i,34:37);
FBAOrd127labels{i} = Convert(i,34:37);
end

countFBAs =0;
for i=1:s2(1)
    if(ReducedFBASolu{i,2}<length(Data.fbas))%if all the FBAs are<128 
          countFBAs= countFBAs+1; %with 'FOR','IA','NA','REV'status =31
            if (ReducedFBASolu{i,3}==UniStates(1,1)&& strcmp(ReducedFBASolu{i,4},'NA'))
                MS_MatrixRow= find(ReducedIDX == ReducedFBASolu{i,1});
                L = fieldnames(Data.reactions{1,ReducedFBASolu{i,1}}.reaction_fluxes);
                Convert = char(L);
                s3= size(Convert); 
                FBAOrdNo127 = cell(s3(1),1);
                for j= 1:s3(1)
                    FBAOrdNo127(j)=cellstr(Convert(j,34:37));
                end
                
                [ExFBA,ExFBAIDX] = setdiff(cellstr(FBAOrd127),cellstr(FBAOrdNo127),'stable'); %ExFBAs-excluded FBAs in the reaction_fluxes(not given) and ExFBAsIDX - index of the excluded fill them with NA states
                temp =ismember((1:length(Data.fbas)),ExFBAIDX);%GSindex-not given state index in the reaction_fluxes
                InFBAIDX = (find(temp==0))'; %InFBAIDX - given FBAs in the reaction_fluxes this can be FOR/REV/IA
               
                NAStates = fieldnames(Data.reactions{1,ReducedFBASolu{i,1}}.state_conservation);
                
                s4 = size(InFBAIDX);
                s5 = size(ExFBAIDX);
                if(strcmp(NAStates(1),'FOR'))
                       for i = 1:s4(1)
                           MS_Matrix(MS_MatrixRow,InFBAIDX(i,1)) = 1;
                       end
                    if(strcmp(NAStates(2),'NA'))
                       for i = 1:s5(1)
                           MS_Matrix(MS_MatrixRow,ExFBAIDX(i,1)) = 0;
                       end
                    end
                end
              
                
                if(strcmp(NAStates(1),'IA'))
                       for i = 1:s4(1)
                           MS_Matrix(MS_MatrixRow,InFBAIDX(i,1)) = 0;
                       end
                    if(strcmp(NAStates(2),'NA'))
                       for i = 1:s5(1)
                           MS_Matrix(MS_MatrixRow,ExFBAIDX(i,1)) = 0;
                       end
                    end
                end
                                
                if(strcmp(NAStates(1),'NA'))
                       for i = 1:s5(1)
                           MS_Matrix(MS_MatrixRow,ExFBAIDX(i,1)) = 0;
                       end
                    if(strcmp(NAStates(2),'REV'))
                       for i = 1:s4(1)
                           MS_Matrix(MS_MatrixRow,InFBAIDX(i,1)) = 1;
                       end  
                    end
                end
           end
   end
end

for i=1:s2(1)
    if(ReducedFBASolu{i,2}<length(Data.fbas))%if all the FBAs are<128 
         if (ReducedFBASolu{i,3}==UniStates(1,2)&& strcmp(ReducedFBASolu{i,4},'NA'))
                MS_MatrixRow= find(ReducedIDX == ReducedFBASolu{i,1});
                L = fieldnames(Data.reactions{1,ReducedFBASolu{i,1}}.reaction_fluxes);
                Convert = char(L);
                s3= size(Convert); 
                FBAOrdNo127 = cell(s3(1),2);
                for j= 1:s3(1)
                    S = strcat('Data.reactions{1,',num2str(ReducedFBASolu{i,1}),'}.reaction_fluxes','.',L{j,1}) ;
                    Q = eval(S);
                    FBAOrdNo127(j,:)= [cellstr(Convert(j,34:37)),Q(1)];
                end
                
                count1=1;
                count2=1;
                count3=1;
                for j= 1:s3(1)
                  if(strcmp(FBAOrdNo127(j,2),'FOR'))
                    FORState(count1,1) = FBAOrdNo127(j,1);
                    count1 = count1+1;
                  elseif(strcmp(FBAOrdNo127(j,2),'IA'))
                    IAState(count2,1) = FBAOrdNo127(j,1);
                    count2 = count2+1;
                  elseif(strcmp(FBAOrdNo127(j,2),'REV')) 
                    REVState(count3,1) = FBAOrdNo127(j,1);
                    count3 = count3+1;
                  end
                end
                
                [ExNAFBA,ExNAFBAIDX] = setdiff(cellstr(FBAOrd127),cellstr(FBAOrdNo127(:,1)),'stable'); %ExFBAs-excluded FBAs in the reaction_fluxes(not given) and ExFBAsIDX - index of the excluded fill them with NA states
%                 
                [ExFORFBA,ExFORFBAIDX] = setdiff(cellstr(FBAOrd127),cellstr(FORState),'stable');
                temp =ismember((1:length(Data.fbas)),ExFORFBAIDX);%GSindex-not given state index in the reaction_fluxes
                InFORFBAIDX = (find(temp==0))';
%                 
%               [ExIAFBA,ExIAFBAIDX] = setdiff(cellstr(FBAOrd127),cellstr(IAState),'stable');
%               temp =ismember((1:length(Data.fbas)),ExIAFBAIDX);%GSindex-not given state index in the reaction_fluxes
%               InIAFBAIDX = (find(temp==0))'
% 
                [ExREVFBA,ExREVFBAIDX] = setdiff(cellstr(FBAOrd127),cellstr(REVState),'stable');
                temp =ismember((1:length(Data.fbas)),ExREVFBAIDX);%GSindex-not given state index in the reaction_fluxes
                InREVFBAIDX = (find(temp==0))';    
                
                %FOR
                 s3 = size(InFORFBAIDX);
                 for j=1:s3(1)
                 MS_Matrix(MS_MatrixRow,InFORFBAIDX(j,1)) = 1;     
                 end       
%                 
                %IA
%                s3 = size(InIAFBAIDX);
%                for j=1:s3(1)
%                MS_Matrix(MS_MatrixRow,InIAFBAIDX(j,1)) = 0;     
%                end    
%                 
                %NA
                 s3 = size(ExNAFBAIDX);
                 for j=1:s3(1)
                 MS_Matrix(MS_MatrixRow,ExNAFBAIDX(j,1)) = 0;     
                 end
%                 
                %REV
                 s3 = size(InREVFBAIDX);
                 for j=1:s3(1)
                 MS_Matrix(MS_MatrixRow,InREVFBAIDX(j,1)) = 1;     
                 end     
         end
    end
end

for i=1:s2(1)
    if(ReducedFBASolu{i,2}<length(Data.fbas))%if all the FBAs are<128 
         if (ReducedFBASolu{i,3}==UniStates(1,1)&& ((strcmp(ReducedFBASolu{i,4},'FOR'))||(strcmp(ReducedFBASolu{i,4},'IA'))||(strcmp(ReducedFBASolu{i,4},'REV'))))%FBAs<128 and the states FOR,IA,REV
         MS_MatrixRow= find(ReducedIDX == ReducedFBASolu{i,1});
         L = fieldnames(Data.reactions{1,ReducedFBASolu{i,1}}.reaction_fluxes);
                Convert = char(L);
                s3= size(Convert); 
                FBAOrdNo127 = cell(s3(1),1);
                for j= 1:s3(1)
                    FBAOrdNo127(j)=cellstr(Convert(j,34:37));
                end
                
                [ExFBA,ExFBAIDX] = setdiff(cellstr(FBAOrd127),cellstr(FBAOrdNo127),'stable'); %ExFBAs-excluded FBAs in the reaction_fluxes(not given) and ExFBAsIDX - index of the excluded fill them with NA states
                temp =ismember((1:length(Data.fbas)),ExFBAIDX);%GSindex-not given state index in the reaction_fluxes
                InFBAIDX = (find(temp==0))'; %InFBAIDX - given FBAs in the reaction_fluxes this can be FOR/REV/IA
               
                NAStates = fieldnames(Data.reactions{1,ReducedFBASolu{i,1}}.state_conservation);
               
                s4 = size(InFBAIDX);
                s5 = size(ExFBAIDX);
                if(strcmp(NAStates(1),'FOR'))
                       for i = 1:s4(1)
                           MS_Matrix(MS_MatrixRow,InFBAIDX(i,1)) = 1;
                       end
                    if(strcmp(NAStates(2),'NA'))
                       for i = 1:s5(1)
                           MS_Matrix(MS_MatrixRow,ExFBAIDX(i,1)) = 0;
                       end
                    end
                end
              
                
                if(strcmp(NAStates(1),'IA'))
                       for i = 1:s4(1)
                           MS_Matrix(MS_MatrixRow,InFBAIDX(i,1)) = 0;
                       end
                    if(strcmp(NAStates(2),'NA'))
                       for i = 1:s5(1)
                           MS_Matrix(MS_MatrixRow,ExFBAIDX(i,1)) = 0;
                       end
                    end
                end
                                
                if(strcmp(NAStates(1),'NA'))
                       for i = 1:s5(1)
                           MS_Matrix(MS_MatrixRow,ExFBAIDX(i,1)) = 0;
                       end
                    if(strcmp(NAStates(2),'REV'))
                       for i = 1:s4(1)
                           MS_Matrix(MS_MatrixRow,InFBAIDX(i,1)) = 1;
                       end  
                    end
                end                      
         end
    end
end

for i=1:s2(1)
    if(ReducedFBASolu{i,2}<length(Data.fbas))%if all the FBAs are<128 
         if (ReducedFBASolu{i,3}==UniStates(1,2)&& ((strcmp(ReducedFBASolu{i,4},'FOR'))||(strcmp(ReducedFBASolu{i,4},'IA'))||(strcmp(ReducedFBASolu{i,4},'REV'))))%FBAs<128 and the states FOR,IA,REV
            MS_MatrixRow= find(ReducedIDX == ReducedFBASolu{i,1});
             L = fieldnames(Data.reactions{1,ReducedFBASolu{i,1}}.reaction_fluxes);
             Convert = char(L);
             s3= size(Convert); 
             FBAOrdNo127 = cell(s3(1),2);
                for j= 1:s3(1)
                    S = strcat('Data.reactions{1,',num2str(ReducedFBASolu{i,1}),'}.reaction_fluxes','.',L{j,1}) ;
                    Q = eval(S);
                    FBAOrdNo127(j,:)= [cellstr(Convert(j,34:37)),Q(1)];
                end
                
                count1=1;
                count2=1;
                count3=1;
                for j= 1:s3(1)
                  if(strcmp(FBAOrdNo127(j,2),'FOR'))
                    FORState(count1,1) = FBAOrdNo127(j,1);
                    count1 = count1+1;
                  elseif(strcmp(FBAOrdNo127(j,2),'IA'))
                    IAState(count2,1) = FBAOrdNo127(j,1);
                    count2 = count2+1;
                  elseif(strcmp(FBAOrdNo127(j,2),'REV')) 
                    REVState(count3,1) = FBAOrdNo127(j,1);
                    count3 = count3+1;
                  end
                end
                
                [ExNAFBA,ExNAFBAIDX] = setdiff(cellstr(FBAOrd127),cellstr(FBAOrdNo127(:,1)),'stable'); %ExFBAs-excluded FBAs in the reaction_fluxes(not given) and ExFBAsIDX - index of the excluded fill them with NA states
                
                [ExFORFBA,ExFORFBAIDX] = setdiff(cellstr(FBAOrd127),cellstr(FORState),'stable');
                temp =ismember((1:length(Data.fbas)),ExFORFBAIDX);%GSindex-not given state index in the reaction_fluxes
                InFORFBAIDX = (find(temp==0))';
                
%               [ExIAFBA,ExIAFBAIDX] = setdiff(cellstr(FBAOrd127),cellstr(IAState),'stable');
%               temp =ismember((1:length(Data.fbas)),ExIAFBAIDX);%GSindex-not given state index in the reaction_fluxes
%               InIAFBAIDX = (find(temp==0))'

                [ExREVFBA,ExREVFBAIDX] = setdiff(cellstr(FBAOrd127),cellstr(REVState),'stable');
                temp =ismember((1:length(Data.fbas)),ExREVFBAIDX);%GSindex-not given state index in the reaction_fluxes
                InREVFBAIDX = (find(temp==0))';    
                
                %FOR
                s3 = size(InFORFBAIDX);
                for j=1:s3(1)
                 MS_Matrix(MS_MatrixRow,InFORFBAIDX(j,1)) = 1;     
                end       
                
                %IA
%                 s3 = size(InIAFBAIDX);
%                 for j=1:s3(1)
%                  MS_Matrix(MS_MatrixRow,InIAFBAIDX(j,1)) = 0;     
%                 end    
                
                %NA
                s3 = size(ExNAFBAIDX);
                for j=1:s3(1)
                 MS_Matrix(MS_MatrixRow,ExNAFBAIDX(j,1)) = 0;     
                end
                
                %REV
                s3 = size(InREVFBAIDX);
                for j=1:s3(1)
                 MS_Matrix(MS_MatrixRow,InREVFBAIDX(j,1)) = 1;     
                end     
                
                
         end
    end
end

for i=1:s2(1)
    ReducedReactionLabel{i,1} = Data.reactions{1,ReducedFBASolu{i,1}}.id;
end

save('MS_Matrix.mat','MS_Matrix')
save('ReducedReactionLabel.mat','ReducedReactionLabel')

msg1 = 'use MS_Matrix.mat, FBAorder.mat, and ReducedReactionLabel.mat to';
msg2 = 'construct MS_2states.txt file using MS_2states_text.m ...';
disp(msg1);
disp(msg2);
clearvars -except Data Matrix AllDifStates FBAOrd127 FBASolu ReducedFBASolu UniStates FBAOrd127labels

