function oscillation_series = make_oscillation(freq,sr,duration)
%oscillation_series = make_oscillation(freq,sr,duration)
%   function oscillation_series = make_oscillation(freq,sr,duration)
%   function to generate oscillating time series.
%   freq: frequency of the oscillation
%   sr : sampling rate, per second
%   duration = total duration of the flickering series in seconds sec

oscillation_series = zeros(duration*sr,1);

ind = 0;
dur_mat = .1:.1:maxtime;
for i = 1:duration*sr
    while ind<duration*sr
        freq = maxfreq*rand;
        dur = dur_mat(randi(length(dur_mat)));
        oscillation_series(ind+1:ind+dur*sr) = sin(2*pi*freq*(1/sr:1/sr:dur));
        ind = ind+dur*sr;
    end
end

  oscillation_series(duration*sr+1:end) = [];
  