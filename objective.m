function D = objective( X )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%% compute predicted signal
global frame_num FreqEchImg sum_signal_art;
S=zeros(frame_num,1);
for j=1:frame_num
   for i=1:length(X)
      S(j) = S(j)+signal(j/FreqEchImg-X(i));
   end    
end
D = S - sum_signal_art;
end

