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
% FreqEch = FreqEchSimu , DureeSimu = DureeAnalysee,NSondeFluo = NbrSondeFluo,NSondeParIntensite=NombreSondeParIntensite
function s=signal(x) 
    [FreqEchS, FreqEchImg, DureeSimu, NSondeFluo, NSondeParIntensite,...
        TaillePreMarq,TailleSeqMarq, TaillePostMarq, VitessePolymerase,frame_num] = deal(Parameters{:});
    Intensity_for_1_Polym = NSondeFluo/NSondeParIntensite;
    T1 = (TaillePreMarq/VitessePolymerase);
    T2 = (TailleSeqMarq/VitessePolymerase);
    T3 = (TaillePostMarq/VitessePolymerase);
    eps = 1e-5;
    signal =0;
    if x > T1 & x < T1 + T2
        s = Intensity_for_1_Polym*(x-T1)/(T2-T1);
    end
    if x > T1 + T2  & x < T1 + T2 + T3
        s = Intensity_for_1_Polym*(1-eps*(x-T1-T2)/(T3-T1-T2));   
    end
end