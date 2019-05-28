function [t_coef] = get_trans_vector(wave_energy)

t_coef = zeros(1,size(wave_energy,2));

for q = 1:size(wave_energy,2)
    t_coef(q) = trans_coef(region_number,potentials,widths,heights,wave_energy(q),wave_amplitude);
end


end