%% Convolution function 
% input: 1 polymerase start position
% output: 1 signal intensity 

% % This function is transforming a Start Polymerase position in signal intensities versus time passed
% % PolyPosition is the polymerase position inside the simulation (startposition in time * Simulation Sampling Frequency)
% % TaillePreMarq/VitessePolymerase duration of the polymerase inside the Premarq
% % NSondeFluo/NSondeParIntensite Intensity for 1 Polymerase (if NSondeParIntensite < NSondeFluo alors 1P -> more than 1 intensity)
% % So the smaller the NSondeParIntensite the smaller the number of Poly needed (the more difficult to fit the exp values)

% the Fluorescence signal intensity
%       T1  T2        T3
%       |   |          |
%       |   /----------|
%       |  /           |
%       | /            |
% ______|/             |_______
% |
% |
% PolyPosition

% default PolyPosition = 0, FreqEchS = 1 , DureeSimu = 1276 ,NSondeFluo = 128,NSondeParIntensite=20
% FreqEch = FreqEchSimu ,NSondeFluo = NbrSondeFluo,NSondeParIntensite=NombreSondeParIntensite
function [n_signal]=getSignal(PolyPosition, Parameters) 
    [FreqEchS, FreqEchImg, DureeAnalysee, NSondeFluo, NSondeParIntensite,...
        TaillePreMarq,TailleSeqMarq, TaillePostMarq, VitessePolymerase,frame_num] = deal(Parameters{:});
    Intensity_for_1_Polym = NSondeFluo/NSondeParIntensite;
    T1 = round((TaillePreMarq/VitessePolymerase)*FreqEchS);
    T2 = round((TailleSeqMarq/VitessePolymerase)*FreqEchS);
    T3 = round((TaillePostMarq/VitessePolymerase)*FreqEchS);
    PolyPosition = round(PolyPosition);
    first = zeros(1,PolyPosition+T1);
    second = Intensity_for_1_Polym/T2:Intensity_for_1_Polym/T2:Intensity_for_1_Polym;
    third = repmat(Intensity_for_1_Polym,1,T3);
    fourth = [first second third];  % the assembled signal
    len_1_sig = T1+T2+T3; % how many "interval" for 1 signal
    fifth = [];
    if(length(fourth)<len_1_sig+round(frame_num*(FreqEchS/FreqEchImg))+1)
        fifth = zeros(1,round(frame_num*(FreqEchS/FreqEchImg))+1);  % the assembled signal is smaller than the desired simulation duration add some 0 to the end
    end 
    sixth = [fourth fifth];
%     sixth = sixth(1:round(frame_num*FreqEchS/FreqEchImg));  % In all cases keep only the signal of the simulated duration length
%     signal = sixth(len_1_sig:round(FreqEchS/FreqEchImg):len_1_sig+frame_num*round(FreqEchS/FreqEchImg));  % return only the values corresponding to the experimental value times points
    ii=1;
    n_signal = zeros(1,frame_num);
    while ii <= frame_num
        index = round((ii-1)*(FreqEchS/FreqEchImg)+len_1_sig); % image start at len_1_sig, because 
        n_signal(ii) = sixth(index);   
        ii= ii+1;
    end
end