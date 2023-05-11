function output = FIR_filter(filter_length, input_signal, impulse_response)

    if length(input_signal) ~= filter_length || length(impulse_response) ~= filter_length
        disp('FIR filter : filter size and signal mismatch')
        output = NaN ;
        return
    end
    
    if size(input_signal) == [filter_length, 1]
        input_signal = transpose(input_signal) ;
    end
    if size(impulse_response) == [1, filter_length]
        impulse_response = transpose(impulse_response) ;
    end
    
    output = input_signal * impulse_response ;
end