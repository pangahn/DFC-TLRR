figure
hold on
t = 1:10;
Yale = [816.2206	363.5551	280.6962	253.6885	233.4542	229.6135	223.2349	225.6845	223.9973	225.9406];
COIL = [803.9583	402.738	360.1468	323.5674	317.7864	313.2421	301.0274	302.4012	302.9717	296.9667];
UCSD = [383.2646	240.5355	190.2112	178.2365	157.8853	147.837	152.7209	142.2906	150.4351	143.933];

plot(t, Yale,'g-s','linewidth',1.5)
plot(t, COIL,'b--o','linewidth',1.5)
plot(t, UCSD,'r:v','linewidth',1.5)

legend({'ExtYale B', 'COIL20', 'UCSD'}, 'FontSize',12);
xlabel('Number of Partition')
ylabel('Running time of DFC-TLRR (s)')