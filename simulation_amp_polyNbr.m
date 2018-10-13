DataFilePath = 'C:/Users/KevinYeung/Desktop/intern Mont/examples_cvs/';
addpath('C:/Users/KevinYeung/Desktop/intern Mont/examples_cvs/');
data_list=struct2cell(dir(fullfile(DataFilePath)));
file_name_list = data_list(1,3:end);
% DataFileName = '20150210_100X_Movie_Calibration_RI510_128_noT_001_visit_2_c2_CALIB_RF.csv';
%% Model and Parameters 
%--------transcription--------
NbrSondeFluo = 128; % 128 sondes per mRNA (ref publi)
TaillePreMarq = 700; % 700 bases (ref publi)
TailleSeqMarq = 5800; % 2900 bases (ref publi)
TaillePostMarq = 1600; % 1600 bases (ref publi)
DureeProcessing = 100; %  en secondes (100s temps moyen ref publi)
% --------polymerases------- 
EspaceInterPolyMin = 30; % en base (40 bases)
Polym_speed = 67; % average speed bases par seconde (Ref publi)
ProbeByIntensitie_nb = NbrSondeFluo; % how many probes are needed to observe 1 intensity
% -------sample freqence------
NombreRandomIter = 1000; %at least 1000 but 10000 is fine for 1200 secondes et FreqEchSimu1
FreqEchImg = (1/3); % 1/3 image par seconde
% FreqEchSimu = 2; % nb de calcul par seconde %% Attention cette donn?e doit ?tre enti?re. -> int ; no float ! ET PROP a la dur?e des calculs
% DureeSimu = 3000; % 1200 secondes Attention cette valeur ne doit pas ?tre inf?rieure a la dur?e exp?rimentale
%CoeffNormalisation = 1
Gap4TrainDef = 2;

exp_data = {};
simu_data = {};

%% =============Do GA for 8 files to find noise parameter==============
% key parameters in this part:  Nbr_simu_DNA, Generations(GA's stop criteria),
%                                
for data_i = 1:8%length(file_name_list) 
    % -------input file----------
    DataFileName = file_name_list{data_i};
    DataExp = table2array(readtable(DataFileName, 'Delimiter','semi'));  
    DataExpSmooth = transp(DataExp(:,2));  % signal intensity in experiment data
    exp_data{end+1} = DataExpSmooth;
    % -------Parameters----------
    DureeSimu = DataExp(end,1)+DataExp(2,1); % (s)
    frame_num = length(DataExp(:,1));
    DureeSignal = round((TaillePreMarq+TailleSeqMarq + TaillePostMarq) / Polym_speed); % (s)
    DureeAnalysee = DureeSignal + DureeSimu - round(TaillePreMarq / Polym_speed); % (s)
    FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed); % how many interval(possible poly start position) in 1s
    num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed));
    Parameters = {FreqEchSimu, FreqEchImg,DureeAnalysee,NbrSondeFluo,ProbeByIntensitie_nb,TaillePreMarq,...
            TailleSeqMarq, TaillePostMarq, Polym_speed,frame_num,num_possible_poly,EspaceInterPolyMin,DureeSimu};
    % ---------estimate poly number---------
    % no necessary 
    eval_poly_num = [];
    GA_fitness_try_poly = @(x) sum((sumSignal(x,Parameters)-DataExpSmooth).^2);
    for num_p = 30:30:900
        fitmin = 1e+10;
        Estimate_polys = zeros(50,num_p);
        for i = 1:50
            poly_posi = randperm(num_possible_poly,num_p); % randomly choose 25 poly position (test1)
            fit_now = GA_fitness_try_poly(poly_posi);
            if fit_now<fitmin
                fitmin = fit_now;
            end
        %     Trans_positions_art = find(Pattern_poly_art(i,:)==1);
        end
        eval_poly_num(end+1) = fitmin;  
    end
    Nbr_poly_estimate = find(eval_poly_num==min(eval_poly_num))*30;

    % -----------initial population for GA----------
    Nbr_simu_DNA = 500; % number of "chromosome" in GA
%     Nbr_poly_estimate = 270; % different for every data
    Pattern_polys = zeros(Nbr_simu_DNA,num_possible_poly);
    for i = 1:Nbr_simu_DNA
        Pattern_polys(i,randperm(num_possible_poly,Nbr_poly_estimate)) = 1; % randomly choose 25 poly position (test1)
    %     Trans_positions_art = find(Pattern_poly_art(i,:)==1);
    end
    
    % ----------GA----------------
    GA_fitness_exp = @(x) sum((sumSignal(find(x==1),Parameters)-DataExpSmooth).^2);
    options = gaoptimset;
    options = gaoptimset(options,'PopulationType', 'bitstring', 'CreationFcn',@gacreationuniform,'MutationFcn',@mutationuniform,...
        'Display','iter','TolFun',1e-8,'Paretofraction',0.35,'Generations',200);

    options = gaoptimset(options,'PopulationSize', Nbr_simu_DNA);
    options = gaoptimset(options,'InitialPopulation', Pattern_polys);
    % x_GA = ga(GA_fitness_art,length(Pattern_poly_art),[],[],[],[],[],[],[],[],options);
    x_GA_exp = ga(GA_fitness_exp,length(Pattern_polys(1,:)),[],[],[],[],[],[],[],[],options);
    simu_data{end+1} = x_GA_exp;
    save(['simu_data(',num2str(data_i),')_200iters_500population'],'x_GA_exp'); 
end

% %% visualize
% s = simu_exp_data;
% for ii = 1:10:30
%    figure(ii)
%    hold on 
%    e_s = s(ii).exp_signal;
%    plot(e_s,'red')
%    plot(sumSignal(s(ii).poly_posi,Parameters),'blue')
%    plot(sumSignal(s(ii+1).poly_posi,Parameters),'green')
% end

%% Add noise
% study the noise from experiment data and simulation data
% each row of signal_simu is the sumSignal of one simu data of
% corresponding exp data (1 file)
signal_simu = sumSignal(find(simu_data{1}==1),Parameters); 
for ii = 2:length(simu_data) 
    signal_simu(ii,:) = sumSignal(find(simu_data{ii}==1),Parameters);
end 
signal_exp = cell2mat(transpose(exp_data)); % each row is one exp file data
diff_simu_exp = signal_simu-signal_exp;
var_window = 5;
varDiff_simuSignal = zeros(2,length(round(min(signal_simu(:))):var_window:round(max(signal_simu(:))))); % store var of 1 window and signal intensity of this window
idx = 1;
for ii = round(min(signal_simu(:))):var_window:round(max(signal_simu(:)))
    df_in_window = diff_simu_exp(find(signal_simu>ii & signal_simu<ii+var_window));
    varDiff_simuSignal(1,idx) = ii+var_window/2;
    varDiff_simuSignal(2,idx) = var(df_in_window);
    idx = idx+1;
end

% -----visualize noise in terms of signal intensity---------
figure(3)
plot(varDiff_simuSignal(1,:),varDiff_simuSignal(2,:),'o')
plot(signal_simu,(signal_simu-signal_exp),'o')

% -------regression----------
% parameters of noise
x = varDiff_simuSignal(1,:);
y = varDiff_simuSignal(2,:);
p = polyfit(x,y,1);
f = polyval(p,x); 
figure(4)
plot(x,y,'o',x,f,'x') 

% --------add noise to simu data-------
mean_sig = mean(diff_simu_exp(:));
var_sig = polyval(p,signal_simu(4,:));
simu_signal_noise = signal_simu(4,:) + mean_sig + (2*(rand(1,length(signal_simu(4,:)))-0.5)).*var_sig;
figure(5)
plot(simu_signal_noise)
hold on 
plot(signal_exp(4,:),'red')

% compare fitness function
sum((simu_signal_noise-exp_data{4}).^2)
sum((signal_simu(4,:)-exp_data{4}).^2)

% save result for noise 'mean' and 'var'
save(['noise'],'mean_sig','var_sig','p')

%% GA again to estimate artificial data (different amplitude)
% try decovolute artificial data for different amplitude value
% -----create artificial data and add noise---------
% parameter can be changed: 
%       Nbr_poly: number of polymerase in artificial data,
%                 to test the robust against different Nbr_poly
%       amplitude_list: test different amplitude of noise
%       also the parameters in GA     
%
% input:
%   Noise is calculated by 'add Noise' part if you did not run "add Noise" you can load file from previous
%   a = load('noise.mat');
%   p = a.p;
%   mean_sig = a.mean_sig;
% output:
%   artificial data stored in one struct <art_data_diff_amp>
%   simulation data stored in a <simu_data_diff_amp_rep> a 1*80 struct
%   each struct is a result of a simulation with 4 fields: polyNbr,amp,trans_posi,trans_posi_GDy 
%   
    

Nbr_poly = 150;
Pattern_poly_art = zeros(1,num_possible_poly);
Pattern_poly_art(randi(num_possible_poly,1,Nbr_poly)) = 1; % randomly choose 25 poly position (test1)
Pattern_poly_art_diff_amp = Pattern_poly_art; % to store artificial poly
position for different amp
Trans_positions_art = find(Pattern_poly_art==1);
sum_signal_art = sumSignal(Trans_positions_art,Parameters);
% add noise
var_sig = polyval(p,sum_signal_art);
amplitude_list = [0:9]*0.3; % control amplitude of noise
% amplitude_list = [0];
% try artificial data for 5 different amplitude value
art_signal_noise = zeros(length(amplitude_list),length(sum_signal_art)); % to store different amplitude signal
for ii = 1:length(amplitude_list)
    art_signal_noise(ii,:) = sum_signal_art + mean_sig + (amplitude_list(ii)*2*(rand(1,length(sum_signal_art))-0.5)).*var_sig;
end

Pattern_poly_simu_diff_amp = zeros(length(amplitude_list),num_possible_poly); % store simulation poly position 
art_signal_noise_diff_amp = art_signal_noise; %  % to store different amplitude signal
Pattern_poly_simu_diff_amp_rep = {};
% repeat GA for 10times 
for   rep_i = 1:10; % iter number for repeating compute same art_data with different noise 0.1-1
    rep_i
    for ii = 1:length(amplitude_list)

        % prepare initial population 
        Nbr_simu_DNA = 500; % number of "chromosome" in GA
        Nbr_poly_estimate = 140; % different for every data£¬guess polyNbr
        Pattern_polys = zeros(Nbr_simu_DNA,num_possible_poly);
        for i = 1:Nbr_simu_DNA
            Pattern_polys(i,randperm(num_possible_poly,Nbr_poly_estimate)) = 1; % randomly choose poly position 
        %     Trans_positions_art = find(Pattern_poly_art(i,:)==1);
        end

        % ----------GA----------------
        GA_fitness_art_noise = @(x) sum((sumSignal(find(x==1),Parameters)-art_signal_noise(ii,:)).^2); % with noise 
        GA_fitness_art = @(x) sum((sumSignal(find(x==1),Parameters)-sum_signal_art).^2); % without noise
        options = gaoptimset;
        options = gaoptimset(options,'PopulationType', 'bitstring', 'CreationFcn',@gacreationuniform,'MutationFcn',@mutationuniform,...
            'Display','iter','TolFun',1e-8,'Paretofraction',0.35,'Generations',100);

        options = gaoptimset(options,'PopulationSize', Nbr_simu_DNA);
        options = gaoptimset(options,'InitialPopulation', Pattern_polys);

        x_GA_art_noise = ga(GA_fitness_art_noise,length(Pattern_polys(1,:)),[],[],[],[],[],[],[],[],options);
        Pattern_poly_simu_diff_amp(ii,:) = x_GA_art_noise;
        save(['simu_art_data_noise_diff_amp_100iters_500population(',num2str(rep_i-1),')'],'Pattern_poly_simu_diff_amp'); 
%         x_GA_art = ga(GA_fitness_art,length(Pattern_polys(1,:)),[],[],[],[],[],[],[],[],options); % simu no noise
        
    end
    Pattern_poly_simu_diff_amp_rep{rep_i} = Pattern_poly_simu_diff_amp; % repeat simu 9 times for 9 diff amp
end
% transform data structure
simu_data_diff_amp_rep = struct;
art_data_diff_amp = struct;
count_s_a = 1;
for amp_i = 1:10
    for rep_i = 1:8
        amp = amplitude_list(amp_i);
        s_p = find(Pattern_poly_simu_diff_amp_rep{rep_i}(amp_i,:));
        s_p_gd = poly_position_diff_amp_GD_rep{rep_i}{amp_i};
        simu_data_diff_amp_rep(count_s_a).polyNbr = 150;
        simu_data_diff_amp_rep(count_s_a).amp = amp;
        simu_data_diff_amp_rep(count_s_a).trans_posi = s_p;
        simu_data_diff_amp_rep(count_s_a).trans_posi_GDy = s_p_gd;
        count_s_a = count_s_a+1;
    end
end
art_data_diff_amp.polyNbr=length(find(Pattern_poly_art_diff_amp==1));
art_data_diff_amp.trans_posi_art=find(Pattern_poly_art_diff_amp==1);




%% GA again to estimate artificial data
% try decovolute artificial data for different polyNbr
% input:
%   Noise is calculated by 'add Noise' part if you did not run "add Noise" you can load file from previous
%   a = load('noise.mat');
%   p = a.p;
%   mean_sig = a.mean_sig;
% output:
%   artificial data stored in <art_data_diff_polyNbr> with 1*10 struct,
%   structs represent different poly_num
%   
%   simulation data stored in a <simu_data_diff_polyNbr_rep> a 1*80 struct
%   each struct is a result of a simulation with 4 fields: polyNbr,amp,trans_posi,trans_posi_GDy 
%   
simu_data_diff_polyNbr_rep = struct;
art_data_diff_polyNbr = struct;
count_simu_polyNbr = 1;
for ii = 1:10
    %create artificial data for different poly nbr
    Nbr_poly = ii*30;
    Pattern_poly_art = zeros(1,num_possible_poly);
    Pattern_poly_art(randperm(num_possible_poly,Nbr_poly)) = 1; % randomly choose 25 poly position (test1)
    Trans_positions_art = find(Pattern_poly_art==1);
    sum_signal_art = sumSignal(Trans_positions_art,Parameters);
    
    % add noise
    var_sig = polyval(p,sum_signal_art);
    % amplitude_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9,1]; % control amplitude of noise
    amplitude = 0.5;
    % try artificial data for different amplitude value
    art_signal_noise_i = sum_signal_art + mean_sig + (amplitude*2*(rand(1,length(sum_signal_art))-0.5)).*var_sig;
    %art_signal_noise_diff_polyNbr(ii,:) = art_signal_noise_i;

    art_data_diff_polyNbr(ii).Amp = 0.5;
    art_data_diff_polyNbr(ii).polyNbr = ii*30;
    art_data_diff_polyNbr(ii).trans_posi_art = Trans_positions_art;
    art_data_diff_polyNbr(ii).signal_noise = art_signal_noise_i;
  
    for rep_i = 1:8
        % prepare initial population 
        Nbr_simu_DNA = 500; % number of "chromosome" in GA
        Nbr_poly_estimate = ii*30-20; % different for every data
        Pattern_polys = zeros(Nbr_simu_DNA,num_possible_poly);
        for i = 1:Nbr_simu_DNA
            Pattern_polys(i,randperm(num_possible_poly,Nbr_poly_estimate)) = 1; % randomly choose poly position 
        %     Trans_positions_art = find(Pattern_poly_art(i,:)==1);
        end

        % ----------GA----------------
        GA_fitness_art_noise = @(x) sum((sumSignal(find(x==1),Parameters)-art_signal_noise_i).^2);
    %     GA_fitness_art = @(x) sum((sumSignal(find(x==1),Parameters)-sum_signal_art).^2);
        options = gaoptimset;
        options = gaoptimset(options,'PopulationType', 'bitstring', 'CreationFcn',@gacreationuniform,'MutationFcn',@mutationuniform,...
            'Display','iter','TolFun',1e-8,'Paretofraction',0.35,'Generations',100);

        options = gaoptimset(options,'PopulationSize', Nbr_simu_DNA);
        options = gaoptimset(options,'InitialPopulation', Pattern_polys);
        % x_GA = ga(GA_fitness_art,length(Pattern_poly_art),[],[],[],[],[],[],[],[],options);
        x_GA_art_noise = ga(GA_fitness_art_noise,length(Pattern_polys(1,:)),[],[],[],[],[],[],[],[],options);
        
        simu_data_diff_polyNbr_rep(count_simu_polyNbr).polyNbr = Nbr_poly;
        simu_data_diff_polyNbr_rep(count_simu_polyNbr).amp = 0.5;
        simu_data_diff_polyNbr_rep(count_simu_polyNbr).trans_posi = find(x_GA_art_noise==1);
        count_simu_polyNbr = count_simu_polyNbr+1;
    end
    save(['simu_data_noise_diff_polyNbr(polyNbr',num2str(ii*30),')_100iters_500population'],'simu_data_diff_polyNbr_rep'); 
    save(['art_data_noise_diff_polyNbr(polyNbr',num2str(ii*30),')_100iters_500population'],'art_data_diff_polyNbr');
end


%% -----visualize GA result------ 
% ==========comparison of different polyNbr ===========
% blue GA data, red artifical experiment data
% figure(1)
% for ii = 1:10
%     subplot(4,3,ii)
%     plot(sumSignal(find(Pattern_poly_art_diff_polyNbr(ii,:)==1),Parameters),'red')
%     hold on 
%     plot(sumSignal(find(Pattern_poly_simu_diff_polyNbr(ii,:)==1),Parameters))
%     plot(art_signal_noise_diff_polyNbr(ii,:),'black')
%     legend('art data','simu data','art data noise')
%     title(['Nbr Poly=',num2str(ii*30),', amp=0.5'])
% end 

% ----------Measure the difference between art and simu data---------
% use Hausdorff distance
% different poly nbr


% ============comparison of different amplitude================
% figure(1)
% for ii = 1:10
%     subplot(5,2,ii)
%     plot(sumSignal(find(Pattern_poly_art_diff_amp==1),Parameters),'red')
%     hold on 
%     plot(sumSignal(find(Pattern_poly_simu_diff_amp(ii,:)==1),Parameters))
%     plot(art_signal_noise_diff_amp(ii,:),'black')
%     legend('art data','simu data','art data noise')
%     title(['Nbr Poly=150, amp=',num2str(ii*0.1)])
% end 

% ----------Measure the difference between art and simu data---------
% use Hausdorff distance
% different amplitude 

%% Gradient descent after GA
% -------------Yang's GD-----------------
% ========different Amplitude===========
% for every transcription start position (x_GA==1), try to shift it (left
% or right), to optimize the result for different amp
poly_position_diff_amp_GD = {}; % store the result after GD
poly_position_diff_amp_GD_rep = cell(1,8); % store the result after GD for 9 repeat experiment
for rep_i = 1:length(Pattern_poly_simu_diff_amp_rep) % depends on how many repeat simu of GA
    rep_i
    for ii = 1:length(amplitude_list)
        Nbr_poly_estimate = length(find(Pattern_poly_simu_diff_amp_rep{rep_i}(ii,:)==1)); % estimate polyNbr by GA 
        sum_signal_art = art_signal_noise_diff_amp(ii,:);
        GD_y_fitness = @(x) sum((sumSignal(x,Parameters)-sum_signal_art).^2); % x: positions of poly
        shift_window = round(num_possible_poly/Nbr_poly_estimate)+10; % one position can move [-s_w,s_w]
        trans_position_y = find(Pattern_poly_simu_diff_amp_rep{rep_i}(ii,:)==1);
        Min_fit = GD_y_fitness(trans_position_y);
        for posi_i = 1: length(trans_position_y)
            new_tr_p = trans_position_y;
            for j = -shift_window:shift_window
                new_tr_p(posi_i) = trans_position_y(posi_i)+j;
                if new_tr_p(posi_i) <= 0
                    continue
                end
                if GD_y_fitness(new_tr_p)<Min_fit
                    trans_position_y = new_tr_p;
                    Min_fit = GD_y_fitness(trans_position_y);
                end
            end
        end
        poly_position_diff_amp_GD{ii} = trans_position_y;
    end
    poly_position_diff_amp_GD_rep{rep_i} = poly_position_diff_amp_GD;
end

% =========Different PolyNbr==============
% for every transcription start position (x_GA==1), try to shift it (left
% or right), to optimize the result for different polyNbr
matlabpool local 2;
poly_position_diff_polyNbr_GD = {}; % store the result after GD
poly_position_diff_polyNbr_GD_rep = cell(1,8); % store the result after GD for 9 repeat experiment
parfor simu_i = 1:length(simu_data_diff_polyNbr_rep) % depends on how many repeat simu of GA
    simu_i
    tic
    s = simu_data_diff_polyNbr_rep(simu_i);
    Nbr_poly_y = s.polyNbr; % estimate polyNbr by GA 
    sum_signal_art = art_data_diff_polyNbr(Nbr_poly_y/30).signal_noise;
    GD_y_fitness = @(x) sum((sumSignal(x,Parameters)-sum_signal_art).^2); % x: positions of poly
    shift_window = round(num_possible_poly/Nbr_poly_y)*2+10; % one position can move [-s_w,s_w]
    trans_position_y = s.trans_posi;
    Min_fit = GD_y_fitness(trans_position_y);
    for posi_i = 1: length(trans_position_y)
        new_tr_p = trans_position_y;
        for j = -shift_window:shift_window
            new_tr_p(posi_i) = trans_position_y(posi_i)+j;
            if new_tr_p(posi_i) <= 0
                continue
            end
            if GD_y_fitness(new_tr_p)<Min_fit
                trans_position_y = new_tr_p;
                Min_fit = GD_y_fitness(trans_position_y);
            end
        end
    end
    simu_data_diff_polyNbr_rep(simu_i).trans_posi_GDy_sw2 = trans_position_y;
    toc
end

matlabpool close

% =========Different shift_window==============
% for every transcription start position (x_GA==1), try to shift it (left
% or right), to optimize the result for different polyNbr
matlabpool local 2;
simu_data_diff_sw = struct;
simu_data_diff_sw = simu_data_diff_polyNbr_rep(33);
simu_data_diff_sw(2:8) = simu_data_diff_polyNbr_rep(34:40);
for rep_i = 1:1 % depends on how many repeat simu of GA
    rep_i
    tic
    Nbr_poly_y = simu_data_diff_sw(rep_i).polyNbr; % estimate polyNbr by GA 
    simu_data_diff_sw(rep_i).trans_posi_diff_sw = zeros(10,length(simu_data_diff_sw(rep_i).trans_posi));
    sum_signal_art = art_data_diff_polyNbr(Nbr_poly_y/30).signal_noise;
    GD_y_fitness = @(x) sum((sumSignal(x,Parameters)-sum_signal_art).^2); % x: positions of poly
    for sw = 1:1
        shift_window = 10*sw; % one position can move [-s_w,s_w]
        trans_position_y = simu_data_diff_sw(rep_i).trans_posi;
        Min_fit = GD_y_fitness(trans_position_y);
        for posi_i = 1: length(trans_position_y)
            new_tr_p = trans_position_y;
            for j = -shift_window:shift_window
                new_tr_p(posi_i) = trans_position_y(posi_i)+j;
                if new_tr_p(posi_i) <= 0
                    continue
                end
                if GD_y_fitness(new_tr_p)<Min_fit
                    trans_position_y = new_tr_p;
                    Min_fit = GD_y_fitness(trans_position_y);
                end
            end
        end
        simu_data_diff_sw(rep_i).trans_posi_diff_sw(sw,:) = trans_position_y;
    end
    simu_data_diff_sw(rep_i).trans_posi_GDy_sw2 = trans_position_y;
    toc
end

matlabpool close


%% Visualize data
% % check signal plot
% figure(10)
% plot(sumSignal(poly_position_diff_amp_GD{10},Parameters))
% hold on 
% plot(art_signal_noise_diff_amp(10,:),'black')
% % check each position
% figure(5)
% plot(poly_position_diff_amp_GD{ii},1,'Marker','o','color','red')
% hold on 
% plot(find(Pattern_poly_art_diff_amp==1),1,'Marker','x','color','blue')

% ============different amplitude==========
% chech Jaccard distance
% check how many position is exactly/approximately the same
jd_dist = @(a,b) 1-length(intersect(a,b))/length(union(a,b)); % compute Jaccard distance
jd_approximate_dist =  @(ref,b) 1-(length(intersect(unique([ref-1,ref-2,ref,ref+1,ref+2]),unique([b-1,b-2,b,b+1,b+2])))/length(union(unique([ref-1,ref-2,ref,ref+1,ref+2]),unique([b-1,b-2,b,b+1,b+2]))));
jaccard_dist_GA = [];
jaccard_dist_GA_rep = zeros(length(Pattern_poly_simu_diff_amp_rep),10); % store Jaccard distance for repeat GA simulation
jaccard_dist_GD = [];
jaccard_dist_GD_rep = zeros(length(Pattern_poly_simu_diff_amp_rep),10);

jaccard_dist_GA_approx = [];
jaccard_dist_GA_rep_approx = zeros(length(Pattern_poly_simu_diff_amp_rep),10); % store Jaccard distance for repeat GA simulation
jaccard_dist_GD_approx = [];
jaccard_dist_GD_rep_approx = zeros(length(Pattern_poly_simu_diff_amp_rep),10);

for rep_i = 1:length(Pattern_poly_simu_diff_amp_rep)
    poly_position_diff_amp_GD = poly_position_diff_amp_GD_rep{rep_i};
    for ii = 1:length(amplitude_list) 
        % jd dist
        art_posi = find(Pattern_poly_art_diff_amp==1);
        simu_posi_GD = poly_position_diff_amp_GD{ii};
        simu_posi_GA = find(Pattern_poly_simu_diff_amp_rep{rep_i}(ii,:)==1);
        jaccard_dist_GD(ii) = jd_dist(art_posi,simu_posi_GD);
        jaccard_dist_GA(ii) = jd_dist(art_posi,simu_posi_GA);
        % how many position is within certain range
        jaccard_dist_GD_approx(ii) = jd_approximate_dist(art_posi,simu_posi_GD);
        jaccard_dist_GA_approx(ii) = jd_approximate_dist(art_posi,simu_posi_GA);
%         same_position_GD(ii) =length( intersect(poly_position_diff_amp_GD{ii},find(Pattern_poly_art_diff_amp==1)));
%         same_position_range_GD(ii) = length(intersect(poly_position_diff_amp_GD{ii},tr_posi_range)); 
    end
%     same_position_range_GD_rep(rep_i,:) = same_position_range_GD;
    jaccard_dist_GA_rep(rep_i,:) = jaccard_dist_GA;
    jaccard_dist_GD_rep(rep_i,:) = jaccard_dist_GD;
    jaccard_dist_GA_rep_approx(rep_i,:) = jaccard_dist_GA_approx;
    jaccard_dist_GD_rep_approx(rep_i,:) = jaccard_dist_GD_approx;
    
    
end

% ============different polyNbr==========
jaccard_dist_GA_polyNbr_rep = zeros(8,10);
jaccard_dist_GD_polyNbr_rep = zeros(8,10);
jaccard_dist_GA_polyNbr_rep_approx= zeros(8,10);
jaccard_dist_GD_polyNbr_rep_approx = zeros(8,10);
%jd_GA_amp=zeros(8,10);
for ii = 1:length(simu_data_diff_polyNbr_rep) 
    s = simu_data_diff_polyNbr_rep(ii);
    polyNbr = s.polyNbr;
    art_posi = art_data_diff_polyNbr([art_data_diff_polyNbr.polyNbr]==polyNbr).trans_posi_art;
    jaccard_dist_GA_polyNbr_rep(ii) = jd_dist(art_posi,s.trans_posi);
    jaccard_dist_GD_polyNbr_rep(ii) = jd_dist(art_posi,s.trans_posi_GDy);
    jaccard_dist_GA_polyNbr_rep_approx(ii) = jd_approximate_dist(art_posi,s.trans_posi);
    jaccard_dist_GD_polyNbr_rep_approx(ii) = jd_approximate_dist(art_posi,s.trans_posi_GDy);
   
%     s_a = simu_data_diff_amp_rep(ii);
%     art_posi_amp = find(Pattern_poly_art_diff_amp==1);
%     jd_GA_amp(ii) = jd_dist(art_posi_amp,s_a.trans_posi);
%     
end

% ============different shift window==========
% verify the effect of different shift window
% polyNbr = 150, amp=0.5, rep = 5, GAiters = 100, GApop=500,
jd_sw = zeros(5,10);
jd_sw_app = zeros(5,10);
same_posi_sw = zeros(5,10);
app_posi_sw = zeros(5,10);
%jd_GA_amp=zeros(8,10);
for rep_i = 1:5
    for sw_i = 1:10
        s = simu_data_diff_sw(rep_i);
        art_posi = art_data_diff_polyNbr(5).trans_posi_art;
        sw_posi = s.trans_posi_diff_sw(sw_i,:);
%         jd_sw(rep_i,sw_i) = jd_dist(art_posi,sw_posi);
%         jd_sw_app(rep_i,sw_i) = jd_approximate_dist(art_posi,sw_posi);
        same_posi_sw(rep_i,sw_i) = length(intersect(art_posi,sw_posi));
        ref = art_posi;
        app_posi_sw(rep_i,sw_i) = length(intersect(unique([ref-1,ref-2,ref,ref+1,ref+2]),sw_posi));
    end
%     s_a = simu_data_diff_amp_rep(ii);
%     art_posi_amp = find(Pattern_poly_art_diff_amp==1);
%     jd_GA_amp(ii) = jd_dist(art_posi_amp,s_a.trans_posi);
%     
end


