SC = [6352.1 6590.15 ; 9039.7 9220.9; 3615.8 4320.3];
X = categorical({'[Altschuler2019]','[Mukherjee2021]','rTT-SPRT'});
X = reordercats(X,{'[Altschuler2019]','[Mukherjee2021]','rTT-SPRT'});
figure()
bar(X,SC)
grid on
legend('\delta=0.1','\delta=0.05')