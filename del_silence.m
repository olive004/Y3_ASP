% ASP CW4 4.5 delete silence from audio files
% threshold should be 2x1 vector
function audio_short = del_silence(audio, window_size, num_sections)
    % default
    thresh = [0.008, 0.008];
    
    m_audio = movmean(abs(audio(:,1)), window_size);
    % find the last value that's too low before the full maximum value
    [a_max, max_i] = max(m_audio);
    % get last
    i1 = max(find( a_max*thresh(1) > m_audio(1:max_i)));
    % get first
    i2 = min(find( a_max*thresh(2) > m_audio(max_i:end)));
    if num_sections>1
        i2 = max_i+i2;
        for section =2:num_sections
            i3 = min(find( a_max*thresh > m_audio(i2:end)));
            i2 = i2+i3;
        end
        i2 = i2-max_i;
    end

    audio_short = audio(i1:i1+i2,1);


end