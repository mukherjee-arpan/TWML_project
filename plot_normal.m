delta = [0.25*10^-1, 0.17*10^-1, 10^-2, 0.5*10^-2, 10^-3];

figure()
plot(delta,T_avg_SE,'ro-')
hold on
plot(delta,T_avg_med,'bx-')
hold on
plot(delta,rSPRT,'k*-')
grid on
xlabel('$\delta$')
ylabel('Average sample complexity')
legend('[Mukherjee, NeurIPS2021]','[Altschuler JMLR2019]','rTT-SPRT')
