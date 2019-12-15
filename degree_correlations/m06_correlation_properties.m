matrix_correlation_prop=zeros(12,12);
matrix_correlation_prop_p=zeros(12,12);

for i=1:12
    temp_p1=properties{i};
    for j=1:12
        
        
        temp_p2=properties{j};
        gene=intersect(temp_p1(:,1),temp_p2(:,1));
        [a, temp1]=ismember(gene,temp_p1(:,1));
        [a, temp2]=ismember(gene,temp_p2(:,1));
        
        [matrix_correlation_prop(i,j),matrix_correlation_prop_p(i,j)]=corr(temp_p1(temp1(temp1>0),2),temp_p2(temp2(temp2>0),2),'type','Spearman');

        
    end
end