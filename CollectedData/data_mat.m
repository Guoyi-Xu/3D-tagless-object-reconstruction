function [R, P, C] = data_mat(data, channel, tagnum, antnum)

R = zeros(tagnum, antnum);
P = zeros(tagnum, antnum);
C = zeros(tagnum, antnum);

for i = 1:tagnum
    for j = 1:antnum
        [~, ~, R(i, j), P(i, j), C(i, j)] = data_group(data, channel, i, j);
    end
end

% Permutation of columns, due to the non-lexicographic ordering of the
% antenna port indices of real antennas.
% The blanks below are reserved for this purpose.

end
