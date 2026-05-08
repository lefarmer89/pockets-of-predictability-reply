function LatexTableFull(Tab,colString,rowString,format,pval,onesided,varargin)

%%%  Input
%    Tab:       a matrix containing the numerical values of the tables
%    colString: is a cell of strings each containing the column label
%    rowString: is a cell of strings each containing the row label 
%    format:    specifies the format
%    varargin:  is an optional argument containing a structure
%               Str.title is the title of the table
%               Str.caption is the caption of the table
%               Str.sidewaystable where the table must be horizontal

% SEE EXAMPLES at the END OF THE FILE

%%%  Output
%    A text to be pasted in a Latex file

if isempty(varargin)==0
    disp('\newpage')
    Str=varargin{1};
    try Title=Str.title; catch me; Title='Put Title Here'; end
    try Caption=Str.caption; catch me; Caption='Put Caption Here'; end
     if isempty(pval)==0
     TabColumnSpacer=['c ' repmat('r @{}l ',1,size(Tab,2)) ' c'];
     else
     TabColumnSpacer=['c ' repmat('c ',1,size(Tab,2)) ' c'];    
     end
      try 
          sideways=Str.sidewaystable; 
          if isequal(sideways,'yes')
          disp('\begin{sidewaystable}')
          elseif isequal(sideways,'no')
          disp('\begin{table}[h!]')    
          end
      catch me; 
          disp('\begin{table}[h!]')
      end
     disp('\centering')
     disp(['\caption{\footnotesize{\textbf{' Title   '.} ' Caption '.}}'])
     disp(['\resizebox{1\textwidth}{!}{'])
     disp(['\begin{tabular}{' TabColumnSpacer '}'])
     disp(['\\'])
end

spaces=5;
Ndecimal=3;
cr=repmat({' \cr'},length(rowString),1);

  if isempty(pval)==0
     temp={};
     for j=2:length(colString) 
      %temp{j-1}=['\multicolumn{2}{c}{' colString{j} '}'];
      temp{j-1}=colString{j};
     end
  disp(sprintf([' %14s & ' repmat(['%' format(1) 's & '],1,size(Tab,2)) ' %4s'],colString{1},temp{:},cr{1})) 
  else
  disp(sprintf([' %14s & ' repmat(['%' format(1) 's & '],1,size(Tab,2)) ' %4s'],colString{1},colString{2:end},cr{1})) 
  end


  if isempty(pval)
      pval=NaN(size(Tab,1),size(Tab,2));
  end
  % Build each row as a single sprintf format string, then format the
  % numeric values once. Per-cell logic decides whether to emit a NaN
  % placeholder ('--') or a stars/daggers-decorated value.
  if ~onesided
      for h=1:size(Tab,1)
          rowVals = {};
          temp=[];
          for j=1:size(Tab,2)
              if isnan(Tab(h,j))
                  temp=[temp '%s & '];
                  rowVals{end+1} = '--';
              elseif pval(h,j)<0.01
                  temp=[temp '$%' format '^{***}$ & '];
                  rowVals{end+1} = Tab(h,j);
              elseif pval(h,j)<0.05
                  temp=[temp '$%' format '^{**}$ & '];
                  rowVals{end+1} = Tab(h,j);
              elseif pval(h,j)<0.1
                  temp=[temp '$%' format '^{*}$ & '];
                  rowVals{end+1} = Tab(h,j);
              else
                  temp=[temp '$%' format '$ & '];
                  rowVals{end+1} = Tab(h,j);
              end
          end
          disp(sprintf([' %14s & ' temp ' %4s'],rowString{h},rowVals{:},cr{h}))
      end
  else
      for h=1:size(Tab,1)
          rowVals = {};
          temp=[];
          for j=1:size(Tab,2)
              v = Tab(h,j);
              if isnan(v)
                  temp=[temp '%s & '];
                  rowVals{end+1} = '--';
              elseif pval(h,j)<0.01 && v >= 0
                  temp=[temp '$%' format '^{***}$ & '];
                  rowVals{end+1} = v;
              elseif pval(h,j)<0.05 && v >= 0
                  temp=[temp '$%' format '^{**}$ & '];
                  rowVals{end+1} = v;
              elseif pval(h,j)<0.1 && v >= 0
                  temp=[temp '$%' format '^{*}$ & '];
                  rowVals{end+1} = v;
              elseif pval(h,j)<0.01 && v < 0
                  temp=[temp '$%' format '^{\\dagger\\dagger\\dagger}$ & '];
                  rowVals{end+1} = v;
              elseif pval(h,j)<0.05 && v < 0
                  temp=[temp '$%' format '^{\\dagger\\dagger}$ & '];
                  rowVals{end+1} = v;
              elseif pval(h,j)<0.1 && v < 0
                  temp=[temp '$%' format '^{\\dagger}$ & '];
                  rowVals{end+1} = v;
              else
                  temp=[temp '$%' format '$ & '];
                  rowVals{end+1} = v;
              end
          end
          disp(sprintf([' %14s & ' temp ' %4s'],rowString{h},rowVals{:},cr{h}))
      end
  end
% else 
%     for h=1:size(Tab,1)
%        w=Tab(h,:);
%          for i=1:length(w)
%          str=num2str(w(i)); 
%          clear tempCell; for ii=1:length(str); tempCell{ii}=str(ii) ;end
%          %toosmall=find(strcmp(tempCell,{'e'}));
%          PointPos=find(strcmp(tempCell,{'.'}));
%          Int{1,i}=char(tempCell(1:PointPos-1))'; 
%          tempCell=[tempCell cellstr(num2str(zeros(Ndecimal,1)))'];
%           
%            if pval(h,i)<0.01
%            Dec{1,i}=[char(tempCell(PointPos:PointPos+1+Ndecimal-1))' '***']; 
%            elseif pval(h,i)<0.05
%            Dec{1,i}=[char(tempCell(PointPos:PointPos+1+Ndecimal-1))' '**']; 
%            elseif pval(h,i)<0.1
%            Dec{1,i}=[char(tempCell(PointPos:PointPos+1+Ndecimal-1))' '*']; 
%            else
%            Dec{1,i}=[char(tempCell(PointPos:PointPos+1+Ndecimal-1))' '']; 
%            end
%          end
% 
%     Tab2=cell(1,length(w)*2);Tab2(1:2:end)=Int; Tab2(2:2:end)=Dec;
%     temp=repmat(['%' num2str(spaces) 's & '],1,length(w)*2);
%     disp(sprintf([' %14s & ' temp ' %4s'],rowString{h},Tab2{1,:},cr{h}))
%     end

if isempty(varargin)==0
    disp('\end{tabular}}')
       try 
          sideways=Str.sidewaystable; 
          if isequal(sideways,'yes')
          disp('\end{sidewaystable}')
          elseif isequal(sideways,'no')
          disp('\end{table}')    
          end
      catch me; 
          disp('\end{table}')
      end
end
 
% Example I
% Tabel=rand(4,2);
% pval=[];
% 
% LatexOptions.title='Example Table';
% LatexOptions.caption='This Table is just and example';
% LatexOptions.sidewaystable='no';
% 
% Row_lab{1}='First';
% Row_lab{2}='Second';
% Row_lab{3}='Third';
% Row_lab{4}='Fourth';
% 
% Col_lab{1}='Casess'; 
% Col_lab{2}='b'; 
% Col_lab{3}='c'; 
% 
% LatexTableFull(Tabel,Col_lab,Row_lab,'9.2f',pval,LatexOptions)  

% 
% THIS IS THE OUTPUT
%
% \newpage
% \begin{table}[h!]
% \centering
% \caption{\footnotesize{\textbf{Example Table.} This Table is just and example.}}
% \resizebox{1\textwidth}{!}{
% \begin{tabular}{c c c  c}
%          Casess &         b &         c &   \cr
%           First &      0.59 &      0.16 &   \cr
%          Second &      0.02 &      0.18 &   \cr
%           Third &      0.43 &      0.42 &   \cr
%          Fourth &      0.31 &      0.09 &   \cr
% \end{tabular}}
% \end{table}