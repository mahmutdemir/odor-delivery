function cparad_temp = recombine_cparad(cont_parad)
n_of_dil = (length(cont_parad)-2)/2;
n_of_fields = length(fieldnames(cont_parad));
inp_names = fieldnames(cont_parad);   % get input channel names


for i = 1:n_of_fields
    cparad_temp(1).(inp_names{i}) = cont_parad(1).(inp_names{i});
    for j = 1:n_of_dil
        cparad_temp(j+1).(inp_names{i}) = cont_parad(2*j).(inp_names{i});
    end
    cparad_temp(n_of_dil+2).(inp_names{i}) = cont_parad(2*n_of_dil+2).(inp_names{i});
%     cparad_temp(n_of_dil+3).(inp_names{i}) = cont_parad(2*n_of_dil+3).(inp_names{i});
end
