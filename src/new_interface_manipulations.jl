function zero_resistance_check(lines::Vector{LineParams},silent::Bool = true)
    res_vec = [Float64(l.r) for l in lines]
    skip_lines = find(res_vec.==0.0)
    if !isempty(skip_lines) && !silent
        info("r=0 check: \tremoving $(length(skip_lines)) lines")
    end
    return skip_lines
end

function zero_length_check(lines::Vector{LineParams},silent::Bool = true)
    len_vec = [Float64(l.len) for l in lines]
    skip_lines = find(len_vec.==0.0)
    if !isempty(zero_len) && !silent
        info("len=0 check: \tremoving $(length(zero_len)) lines")
    end
    return skip_lines
end

function shift_factor_check(shift_factors,silent::Bool = true; tiny=1e-8)
    skip_lines = find(1 - [maxabs(shift_factors[i,:]) > small for i in 1:size(shift_factors,1)])
    if !isempty(singular_line_idx) && !silent
        info("shift factor check: \tremoving $(length(singular_line_idx)) lines")
    end
    return skip_lines
end
