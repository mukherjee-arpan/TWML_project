mu = [5,4.5,1,1,1];
e = 0.1;
max_iter = 100000;
K = length(mu);
delta = [0.25*10^-1, 0.17*10^-1, 10^-2, 0.5*10^-2, 10^-3];
avg = 200;
T_avg_med = zeros(1,length(delta));
N_avg_med = zeros(length(delta),K);
success_med = zeros(1,length(delta));
for delta_iter = 1:length(delta)
    for iter  = 1:avg
        iter
        t = 1;
        S = [1,2,3,4,5];
        mu_hat = zeros(1,K);
        Na = zeros(1,K);
        alpha = sqrt(2*log(K*pi^2/(6*delta(delta_iter))));
        data = zeros(K,max_iter);
        while length(S)>1
            for i=1:length(S)
                D = binornd(1,e);
                if D == 1
                    x = 0.3 + rand();
                else
                    x = randn()+ mu(S(i));
                end
                if t == 1
                    mu_hat(S(i)) = x;
                    data(S(i),Na(S(i))+1) = x;
                    Na(S(i)) = Na(S(i)) + 1;
                    t = t+1;
                else
                    data(S(i),Na(S(i))+1) = x;
                    Na(S(i)) = Na(S(i)) + 1;
                    mu_hat(S(i)) = median(data(S(i),1:Na(S(i))));
                    t = t+1;
                end
            end
            elem = [];
            for i=1:length(S)
                if mu_hat(i) <= max(mu_hat) - 2*alpha 
                    elem = [S(i), elem];
                end
            end
            S = setdiff(S,elem);
            alpha = sqrt((2/t)*log(K*t^2*pi^2/(6*delta(delta_iter))));
        end
        success_med(delta_iter) = success_med(delta_iter) + (S(1)==1);
        T_avg_med(delta_iter) = T_avg_med(delta_iter) + t;
        N_avg_med(delta_iter,:) = N_avg_med(delta_iter,:) + Na;
    end
end
success_med(delta_iter) = success_med(delta_iter)/avg;
T_avg_med(delta_iter) = T_avg_med(delta_iter)/avg;
N_avg_med = N_avg_med/avg;
alloc_med = N_avg_med/sum(N_avg_med);