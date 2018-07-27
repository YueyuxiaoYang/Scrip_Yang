%% Sum signals
% ------core-----
% input: transcription start positions(Trans_position), Parameters
% output: sum of signals
% call function: getSignal()
% Parameters = {FreqEchSimu, FreqEchImg,DureeSimu,NbrSondeFluo,...
%            ProbeByIntensitie_nb,TaillePreMarq,TailleSeqMarq, TaillePostMarq, Polym_speed};

function [Sum_signals] = sumSignal(Trans_positions,Parameters)
    [FreqEchSimu, FreqEchImg,DureeSimu,NbrSondeFluo,ProbeByIntensitie_nb,...
            TaillePreMarq,TailleSeqMarq, TaillePostMarq, Polym_speed,frame_num,num_possible_poly,EspaceInterPolyMin,DureeSimu] = deal(Parameters{:});
    l_signal = length(getSignal(1, Parameters));
    Sum_signals = zeros(1,l_signal);
    signals = [];
    for posi_i = 1:length(Trans_positions)
        signal = getSignal(Trans_positions(posi_i), Parameters);
%         signals(posi_i,1:length(signal)) = signal;
        Sum_signals(end-length(signal)+1:end) = Sum_signals(end-length(signal)+1:end) +signal;
    end         
%     % -----plot signals------
%     for ii = 1:length(signals(:,1))
%         plot(signals(ii,:));
%         hold on;
%     end

end 