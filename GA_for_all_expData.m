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

%% Do GA for 8 files to find noise parameter
for data_i = 1:8%length(file_name_list) 
    % -------input file----------
    DataFileName = file_name_list{data_i}
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
    Parameters = {FreqEchSimu, FreqEchImg,DureeAnalysee,NbrSondeFluo,...
            ProbeByIntensitie_nb,TaillePreMarq,TailleSeqMarq, TaillePostMarq, Polym_speed,frame_num};
    % ---------estimate poly number---------
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
    Nbr_poly_estimate = find(eval_poly_num==min(eval_poly_num))*30

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

%% visualize
for ii = 5:8
   figure(ii)
   plot(sumSignal(find(simu_data{ii}==1),Parameters))
   hold on 
   plot(exp_data{ii},'red')
end

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

sum((simu_signal_noise-exp_data{4}).^2)
sum((signal_simu(4,:)-exp_data{4}).^2)



%% GA again to estimate artificial data (different amplitude)
% try decovolute artificial data for different amplitude value
% -----create artificial data and add noise---------
Nbr_poly = 150;
Pattern_poly_art = zeros(1,num_possible_poly);
Pattern_poly_art(randi(num_possible_poly,1,Nbr_poly)) = 1; % randomly choose 25 poly position (test1)
Pattern_poly_art_diff_amp = Pattern_poly_art; % to store artificial poly position for different amp
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
        Nbr_poly_estimate = 140; % different for every data��guess polyNbr
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

%% GA again to estimate artificial data
% try decovolute artificial data for different polyNbr
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
figure(1)
for ii = 1:10
    subplot(4,3,ii)
    plot(sumSignal(find(Pattern_poly_art_diff_polyNbr(ii,:)==1),Parameters),'red')
    hold on 
    plot(sumSignal(find(Pattern_poly_simu_diff_polyNbr(ii,:)==1),Parameters))
    plot(art_signal_noise_diff_polyNbr(ii,:),'black')
    legend('art data','simu data','art data noise')
    title(['Nbr Poly=',num2str(ii*30),', amp=0.5'])
end 


Pattern_poly_art_no_noise = x_GA_art;
Pattern_poly_list_diff_amp = Pattern_poly_list; % store data for different amp value, load data from file
distance_pattern_poly = @(x) sum((x(1:length(Trans_positions_art))-Trans_positions_art).^2);
figure(1)
plot(find(Pattern_poly_list(1,:)==1),1,'Marker','o','color','red')
hold on 
plot(Trans_positions_art,1,'Marker','x','color','blue')

% ----------Measure the difference between art and simu data---------
% use Hausdorff distance
% different poly nbr


% ============comparison of different amplitude================
figure(1)
for ii = 1:10
    subplot(5,2,ii)
    plot(sumSignal(find(Pattern_poly_art_diff_amp==1),Parameters),'red')
    hold on 
    plot(sumSignal(find(Pattern_poly_simu_diff_amp(ii,:)==1),Parameters))
    plot(art_signal_noise_diff_amp(ii,:),'black')
    legend('art data','simu data','art data noise')
    title(['Nbr Poly=150, amp=',num2str(ii*0.1)])
end 

% ----------Measure the difference between art and simu data---------
% use Hausdorff distance
% different amplitude 



% check how many position is exactly the same
jd_dist = @(a,b) length(intersect(a,b))/length(union(a,b)); % compute Jaccard distance
jaccard_dist_GA = [];
jaccard_dist_GA_rep = zeros(length(Pattern_poly_simu_diff_amp_rep),10);
same_position_range_GA = [];
same_position_range_GA_rep=zeros(length(Pattern_poly_simu_diff_amp_rep),10); % 9 GA result for same art_data
tr_posi= find(Pattern_poly_art_diff_amp==1);
tr_posi_range = [tr_posi,tr_posi+1,tr_posi+2,tr_posi-1,tr_posi-2]; % permit_range = 2;
for rep_i = 2:10
    for ii = 1:10
        jaccard_dist_GA(ii) =length( intersect(find(Pattern_poly_simu_diff_amp_rep{rep_i}(ii,:)==1),find(Pattern_poly_art_diff_amp==1)));
        same_position_range_GA(ii) = length(intersect(find(Pattern_poly_simu_diff_amp_rep{rep_i}(ii,:)==1),tr_posi_range));
    end
    same_position_range_GA_rep(rep_i-1,:)=same_position_range_GA;
end
boxplot(same_position_range_GA_rep,'color','b')
plot(mean(same_position_range_GA_rep),'o')
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
    shift_window = round(num_possible_poly/Nbr_poly_y)+10; % one position can move [-s_w,s_w]
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
    simu_data_diff_polyNbr_rep(simu_i).trans_posi_GDy = trans_position_y;
    toc
end

matlabpool close


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

% chech hausdorff distance
% check how many position is exactly/approximately the same
jd_dist = @(a,b) length(intersect(a,b))/length(union(a,b)); % compute Jaccard distance
jd_approximate_dist =  @(ref,b) length(intersect(unique([ref-1,ref-2,ref,ref+1,ref+2]),b))/length(union(ref,b));
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

% visualize 9 repeat result 
figure(3)
boxplot(jaccard_dist_GA_rep,'color','r')
hold on
boxplot(jaccard_dist_GD_rep,'color','b')
plot(mean(jaccard_dist_GA_rep),'marker','o','color','r')
plot(mean(jaccard_dist_GD_rep),'marker','o','color','b')
legend('Jaccard distance after GA','Jaccard distance after GD')
title('Compare Jaccard distance of 10 amplitudes(0-2.7) with 8 repeated simulation(same art_data)')
xlabel('Amplitude of noise(x0.3)')
ylabel('Jaccard distance')

figure(4)
boxplot(jaccard_dist_GA_rep_approx,'color','r')
hold on
boxplot(jaccard_dist_GD_rep_approx,'color','b')
plot(mean(jaccard_dist_GA_rep_approx),'marker','o','color','r')
plot(mean(jaccard_dist_GD_rep_approx),'marker','o','color','b')
legend('Jaccard distance after GA','Jaccard distance after GD')
title('Compare approximate Jaccard distance of 10 amplitudes(0-2.7) with 8 repeated simulation(same art_data)')
xlabel('Amplitude of noise(x0.3)')
ylabel('Approximate Jaccard distance')

figure(4)
hold on 
boxplot(same_position_range_GD_rep,'color','b')
boxplot(same_position_range_GA_rep,'color','r')
plot(mean(same_position_range_GA_rep),'marker','o','color','r')
plot(mean(same_position_range_GD_rep),'marker','o','color','b')
ylim([20,100])
legend('approximately same position nbr after GA','approximately same position nbr after GD')
title('Compare approxiamtely same position nbr of 9 different amplitude with 9 repeated simulation(150 poly)')
xlabel('Amplitude of noise(x0.1)')
ylabel('approximately same position nbr')

figure(3)
plot(hd_dist_list_GD_rep(2,:),'marker','o','color','blue')
hold on
plot(hd_dist_list,'marker','o','color','red')
plot(same_position_range_GD,'marker','x','color','black')
plot(same_position_range_GA,'marker','x','color','green')
legend('hd distance after GD','hd distance after GA','approximately same position nbr after GD','approximately same position nbr after GA')



% check how many position is exactly the same
% check polyNbr 

length(poly_position_diff_amp_GD{ii})
same_position = intersect(poly_position_diff_amp_GD{ii},find(Pattern_poly_art_diff_amp==1));
% check how many position is no more than permit_range
tr_posi = find(Pattern_poly_art_diff_amp==1);
tr_posi_range = [tr_posi,tr_posi+1,tr_posi+2,tr_posi-1,tr_posi-2]; % permit_range = 2;
same_position_range2 = intersect(poly_position_diff_amp_GD{ii},tr_posi_range);


% try parallel computing

c = 1:10000000;
a = ones(10000000,1);
tic
for i = 1:length(c)
    a(i)= a(i)+ c(i);
end
toc

matlabpool local 2;
c = 1:10000000;
a = ones(10000000,1);
tic
parfor i = 1:length(c)
    a(i)= a(i)+ c(i);
end
toc

matlabpool close



