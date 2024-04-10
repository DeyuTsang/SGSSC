function res = partition(v,k)
% PARTITION Partitions a vector v into k contiguous pieces of nearly
% equal length.
num_per_part = floor(length(v)/k);
num_long_parts = length(v) - k*num_per_part;
res = {};
part_start = 1;
for ii = 1:k
   if ii <= num_long_parts
      part_len = num_per_part+1;
   else
      part_len = num_per_part;
   end
   res{ii} = v(part_start:(part_start+part_len-1));
   part_start = part_start+part_len;
end
end
