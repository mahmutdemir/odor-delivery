function binary_series = make_pseudo_binary_series(corr_length,sr,duration)
%make_pseudo_binary_seris(corr_length,sr,duration)
%   function binary_series = make_pseudo_binary_seris(corr_length,sr,duration)
%   function to generate pseudo random binary series for flickering odor
%   stimulus valve control
%   corr_length = minimum length of the time that valve is open or closed in ms
%   sr : sampling rate, per second
%   duration = total duration of the flickering series in seconds sec

binary_series = zeros(duration*sr,1);

ind = 0;
for i = 1:duration*sr
    while ind<duration*sr
        if rand > .5
            binary_series(ind+1:ind+corr_length/1000*sr) = 1;
        else
            binary_series(ind+1:ind+corr_length/1000*sr) = 0;
        end
        ind = ind+corr_length/1000*sr;
    end
end

  binary_series(duration*sr+1:end) = [];
  