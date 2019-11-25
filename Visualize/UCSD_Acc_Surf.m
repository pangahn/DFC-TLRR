X = [0.1 0.5:0.5:15];
Y = [0 0.1:0.1:0.9 0.99];

Z = [
0.8181	0.7254	0.7199	0.721	0.7221	0.7232	0.7232	0.7232	0.7232	0.7243	0.7243	0.7243	0.7266	0.7266	0.7266	0.7266	0.7266	0.7236	0.7266	0.7277	0.7277	0.7277	0.7277	0.808	0.8025	0.8025	0.8025	0.8025	0.7935	0.8036	0.8025
0.8203	0.8225	0.7243	0.7221	0.7221	0.7232	0.7232	0.7232	0.7243	0.7232	0.7232	0.7254	0.7254	0.7734	0.7779	0.7779	0.7812	0.78549	0.7824	0.7846	0.7891	0.7891	0.808	0.8069	0.808	0.8069	0.808	0.808	0.8069	0.808	0.8069
0.8203	0.8225	0.7243	0.7221	0.7232	0.7232	0.7232	0.7232	0.7243	0.7243	0.7243	0.7243	0.7254	0.7734	0.7734	0.7779	0.779	0.7827	0.7824	0.7801	0.7835	0.7835	0.7868	0.7879	0.7935	0.7958	0.798	0.7969	0.798	0.798	0.7991
0.8214	0.8237	0.7266	0.721	0.7232	0.7232	0.7232	0.7232	0.7232	0.7243	0.7232	0.7254	0.7254	0.7746	0.7746	0.7746	0.7746	0.7782	0.7779	0.7812	0.7812	0.7824	0.7812	0.7835	0.7857	0.7879	0.7868	0.7891	0.7879	0.7958	0.7969
0.8248	0.8281	0.7232	0.7221	0.7221	0.7232	0.7232	0.7232	0.7232	0.7243	0.7254	0.7254	0.7746	0.7746	0.7746	0.7746	0.7757	0.7736	0.7757	0.7757	0.7757	0.7757	0.7757	0.7757	0.779	0.7812	0.779	0.7779	0.7779	0.7846	0.7857
0.827	0.8304	0.7232	0.7221	0.7221	0.7232	0.7232	0.7243	0.7243	0.7243	0.7254	0.7266	0.7746	0.7746	0.7746	0.7757	0.8013	0.8036	0.8058	0.8058	0.8058	0.8058	0.8058	0.8058	0.808	0.808	0.808	0.808	0.7991	0.7958	0.7913
0.8292	0.8337	0.8348	0.7667	0.7679	0.7679	0.769	0.769	0.769	0.7701	0.7712	0.7757	0.7779	0.7857	0.7924	0.7924	0.7991	0.8058	0.8036	0.8036	0.8058	0.8092	0.8839	0.885	0.875	0.885	0.8839	0.8817	0.8817	0.8817	0.8817
0.8571	0.8504	0.8482	0.8549	0.856	0.8717	0.8806	0.8795	0.8783	0.8795	0.8806	0.885	0.8884	0.8929	0.8973	0.9085	0.9118	0.8795	0.9118	0.9118	0.9118	0.9107	0.9051	0.8839	0.8839	0.8783	0.8839	0.8783	0.8783	0.8795	0.8783
0.9152	0.7768	0.7846	0.7969	0.8225	0.8103	0.8304	0.8962	0.8961	0.8917	0.8951	0.8951	0.875	0.8962	0.8951	0.8951	0.8951	0.8928	0.8929	0.894	0.894	0.8917	0.8929	0.8906	0.8906	0.8906	0.8906	0.9051	0.8906	0.9051	0.8906
0.8125	0.8158	0.8058	0.8147	0.7991	0.798	0.8002	0.8058	0.8125	0.8136	0.8136	0.8147	0.8125	0.8147	0.8125	0.8125	0.8125	0.7545	0.8103	0.8092	0.7891	0.7879	0.7824	0.7824	0.7835	0.7835	0.7891	0.7879	0.7913	0.7913	0.7913
0.7254	0.7264	0.769	0.7254	0.7623	0.7384	0.615	0.6138	0.6083	0.6462	0.6239	0.6507	0.6283	0.5201	0.5212	0.5558	0.5201	0.635	0.5636	0.529	0.5413	0.5156	0.5201	0.5658	0.5681	0.5681	0.558	0.5714	0.5714	0.5703	0.5569
];


xmin=min(X);xmax=max(X);
ymin=min(Y);ymax=max(Y);

[lx,ly,lz] = griddata(X, Y, Z, linspace(xmin,xmax)',linspace(ymin,ymax),'cubic');

figure
subplot(1, 2, 1)
surf(lx,ly,lz)
xlabel('\lambda')
xticks(1:15)
ylabel('S')
yticks(0:0.2:1)
zlabel('Accuracy')

subplot(1, 2, 2)
pcolor(X,Y,Z)
xlabel('\lambda')
xticks(1:15)
ylabel('S')
yticks(0:0.2:1)
shading interp
zmax=max(max(Z));
zmin=min(min(Z));
caxis([zmin,zmax])
colorbar
hold on
C=contour(X,Y,Z,'k:');
% clabel(C)
hold off