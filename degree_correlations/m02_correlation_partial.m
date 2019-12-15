matrix_partial_correlation=zeros(9,11);
matrix_partial_correlation_p=zeros(9,11);
control=10;
control_properties=properties{control};
for i=1:9
    temp_degree=degree{i};
    for j=1:11
        if j>control-1
            j_properties=j+1;
        else
            j_properties=j;
        end
        %j_properties
        temp_properties=properties{j_properties};
        gene=intersect(intersect(temp_degree(:,1),temp_properties(:,1)),control_properties(:,1));
        
        [a, temp1]=ismember(gene,temp_degree(:,1));
        [a, temp2]=ismember(gene,temp_properties(:,1));
        [a, temp3]=ismember(gene,control_properties(:,1));
        [matrix_partial_correlation(i,j),matrix_partial_correlation_p(i,j)]=partialcorr(temp_properties(temp2(temp2>0),2),temp_degree(temp1(temp1>0),2),control_properties(temp3(temp3>0),2),'type','Spearman');
        %[matrix_partial_correlation(i,j),matrix_partial_correlation_p(i,j)]=partialcorr(temp_properties(temp2(temp2>0),2),temp_degree(temp1(temp1>0),2),control_properties(temp3(temp3>0),2));

    end
end
