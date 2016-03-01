function data_temp = recombine_data(data)
if mod(length(data),2)==0
    n_of_dil = (length(data)-2)/2;
    n_of_fields = length(fieldnames(data));
    inp_names = fieldnames(data);   % get input channel names

    [n_of_trial, dlen] =  size(data((length(data)-2)).(inp_names{1}));



    for i = 1:n_of_fields
        data_temp(1).(inp_names{i}) = data(1).(inp_names{i});
        for j = 1:n_of_dil
            data_temp(j+1).(inp_names{i})(1:n_of_trial,1:dlen) = data(2*j).(inp_names{i});
            data_temp(j+1).(inp_names{i})(n_of_trial+1,1:dlen) = data(2*j+1).(inp_names{i});
        end
         data_temp(n_of_dil+2).(inp_names{i}) = data(2*n_of_dil+2).(inp_names{i});
    %     data_temp(n_of_dil+3).(inp_names{i}) = data(2*n_of_dil+3).(inp_names{i});
    end

else
    n_of_dil = (length(data)-1)/2;
    n_of_fields = length(fieldnames(data));
    inp_names = fieldnames(data);   % get input channel names

    [n_of_trial, dlen] =  size(data((length(data)-1)).(inp_names{1}));

    for i = 1:n_of_fields
        data_temp(1).(inp_names{i}) = data(1).(inp_names{i});
        for j = 1:n_of_dil
            data_temp(j+1).(inp_names{i})(1:n_of_trial,1:dlen) = data(2*j).(inp_names{i});
            data_temp(j+1).(inp_names{i})(n_of_trial+1,1:dlen) = data(2*j+1).(inp_names{i});
        end
    %     data_temp(n_of_dil+2).(inp_names{i}) = data(2*n_of_dil+2).(inp_names{i});
    %     data_temp(n_of_dil+3).(inp_names{i}) = data(2*n_of_dil+3).(inp_names{i});
    end
end
