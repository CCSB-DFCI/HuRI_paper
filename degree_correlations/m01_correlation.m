matrix_correlation=zeros(9,12);
matrix_correlation_p=zeros(9,12);

for i=1:9
    temp_degree=degree{i};
    for j=1:12
        
        temp_properties=properties{j};
        gene=intersect(temp_degree(:,1),temp_properties(:,1));
        [a, temp1]=ismember(gene,temp_degree(:,1));
        [a, temp2]=ismember(gene,temp_properties(:,1));
        
        [matrix_correlation(i,j),matrix_correlation_p(i,j)]=corr(temp_degree(temp1(temp1>0),2),temp_properties(temp2(temp2>0),2),'type','Spearman');
        %[matrix_correlation(i,j),matrix_correlation_p(i,j)]=corr(temp_degree(temp1(temp1>0),2),temp_properties(temp2(temp2>0),2));
        
    end
end