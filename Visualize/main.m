Y = [0 0.1:0.1:0.9 0.99];

% COIL
Z_COIL = [
0.7049	0.6549	0.6486	0.6875	0.6465	0.6479	0.6785	0.6965	0.6618	0.6889	0.6736	0.641	0.6681	0.6424	0.6403	0.6389	0.6382	0.6306	0.6382	0.6382	0.6382	0.6354	0.6375	0.6361	0.6375	0.6375	0.6389	0.6403	0.6257	0.6403	0.6778	0.6403	0.6743	0.641	0.6778	0.6653	0.6188	0.6556	0.6667	0.6375
0.6639	0.6771	0.6681	0.6667	0.6882	0.6764	0.6458	0.6583	0.6701	0.6833	0.6965	0.6868	0.6819	0.6729	0.6722	0.6722	0.6722	0.675	0.6736	0.675	0.6743	0.6819	0.6715	0.6625	0.6729	0.6944	0.6729	0.6674	0.6632	0.6569	0.6604	0.6556	0.6576	0.6403	0.6618	0.6514	0.6604	0.6507	0.6618	0.6764
0.7181	0.675	0.6729	0.6757	0.6833	0.7007	0.6667	0.6701	0.6681	0.6743	0.6653	0.6472	0.6764	0.6806	0.6417	0.6722	0.6632	0.6771	0.6604	0.675	0.6562	0.6556	0.6708	0.6514	0.6722	0.6826	0.6722	0.6514	0.6681	0.6826	0.6521	0.6701	0.6778	0.6264	0.6681	0.6667	0.675	0.6667	0.6729	0.6653
0.6972	0.6687	0.6576	0.6528	0.6611	0.6757	0.6778	0.6813	0.6687	0.6694	0.6708	0.6757	0.6729	0.6639	0.6549	0.6674	0.691	0.6889	0.6549	0.6729	0.6826	0.6708	0.6549	0.6771	0.6549	0.6799	0.6667	0.6799	0.6674	0.6674	0.6743	0.6847	0.6778	0.6924	0.6819	0.6625	0.659	0.6479	0.6507	0.65
0.5972	0.6576	0.6569	0.6625	0.6493	0.6444	0.6535	0.6389	0.6549	0.6431	0.6882	0.6792	0.6875	0.675	0.6958	0.6549	0.6694	0.6646	0.6514	0.6778	0.7014	0.6514	0.6701	0.6965	0.6361	0.6917	0.6771	0.6951	0.6799	0.6757	0.6757	0.6993	0.6854	0.6757	0.6653	0.675	0.675	0.675	0.6625	0.6757
0.5799	0.5431	0.616	0.6507	0.6368	0.6625	0.6194	0.6556	0.6757	0.6528	0.6937	0.6729	0.6556	0.6465	0.6403	0.6438	0.6556	0.641	0.6389	0.6535	0.6507	0.6618	0.6535	0.6361	0.659	0.6583	0.6382	0.666	0.6389	0.6368	0.6507	0.6549	0.6292	0.6451	0.6514	0.6139	0.6167	0.6451	0.6458	0.6097
0.5611	0.5278	0.5396	0.5819	0.6438	0.6708	0.6569	0.6521	0.6306	0.6424	0.6757	0.6639	0.6882	0.6507	0.6694	0.6868	0.666	0.6549	0.6319	0.6438	0.6312	0.6562	0.6097	0.6382	0.6285	0.6056	0.6264	0.609	0.6278	0.6368	0.6361	0.6021	0.6174	0.6382	0.6097	0.6368	0.6438	0.6062	0.6153	0.641
0.5431	0.5264	0.5326	0.5472	0.5729	0.575	0.6194	0.6194	0.6132	0.6167	0.6125	0.6167	0.5896	0.6056	0.616	0.6257	0.6243	0.5882	0.6021	0.5917	0.6049	0.6118	0.5868	0.6139	0.5917	0.5764	0.5757	0.6146	0.5785	0.5833	0.609	0.5889	0.5868	0.5687	0.5924	0.559	0.6049	0.5847	0.591	0.5951
0.5444	0.5069	0.5118	0.5521	0.5222	0.5292	0.5618	0.5681	0.616	0.6083	0.6229	0.5889	0.5889	0.6264	0.6167	0.5778	0.6271	0.6097	0.5986	0.6118	0.5979	0.6118	0.6056	0.6028	0.6104	0.5958	0.6062	0.6236	0.6049	0.609	0.5715	0.6042	0.6312	0.6028	0.5951	0.6229	0.5965	0.6056	0.5979	0.5944
0.6097	0.5451	0.5944	0.6	0.6236	0.6306	0.6333	0.6507	0.6792	0.6521	0.6736	0.666	0.6708	0.6736	0.675	0.6549	0.6521	0.6458	0.6618	0.6549	0.6403	0.6542	0.6424	0.6604	0.666	0.6375	0.6431	0.6396	0.6708	0.6368	0.659	0.6604	0.6542	0.6424	0.6542	0.6382	0.6646	0.6361	0.6715	0.6431
0.7389	0.7319	0.7014	0.7264	0.7375	0.7472	0.7389	0.7472	0.7861	0.784	0.7868	0.7703	0.7721	0.7771	0.7743	0.7819	0.7799	0.7806	0.7826	0.7715	0.7743	0.7729	0.7806	0.7729	0.7708	0.7743	0.7681	0.7764	0.7743	0.7729	0.775	0.7743	0.775	0.7736	0.766	0.775	0.7757	0.7764	0.7451	0.7438
];
X_COIL = [0.1	0.5	1	1.5	2	2.5	3	3.5	4	4.5	5	5.5	6	6.5	7	8	8.5	9	9.5	10	10.5	11	11.5	12	12.5	13	13.5	14	14.5	15	15.5	16	16.5	17	17.5	18	18.5	19	19.5	20];

xmin=min(X_COIL);xmax=max(X_COIL);
ymin=min(Y);ymax=max(Y);

[lx,ly,lz] = griddata(X_COIL, Y, Z_COIL, linspace(xmin,xmax)',linspace(ymin,ymax),'cubic');

figure

subplot(1, 3, 1)

surf(lx,ly,lz)
title('Acc(\it{s}, \lambda) on COIL20')
xlabel('\lambda')
xticks(1:2:20)
ylabel('\it{s}')
yticks(0:0.2:1)
zlabel('Accuracy')


% Yale B
X_Yale = [0.1 0.5:0.5:8];

Z_Yale = [
0.9414	0.9509	0.9556	0.9643	0.9651	0.9675	0.9865	0.9643	0.9913	0.9913	0.9913	0.9913	0.9612	0.9913	0.9913	0.9913	0.9913
0.9493	0.9881	0.9905	0.9905	0.9897	0.9913	0.9913	0.9913	0.9913	0.9905	0.9723	0.9913	0.9913	0.9913	0.9762	0.9723	0.9731
0.8938	0.9715	0.9659	0.9683	0.9913	0.9897	0.9905	0.9699	0.912	0.9715	0.9136	0.9715	0.912	0.9723	0.9128	0.912	0.912
0.8947	0.9057	0.9097	0.9089	0.912	0.9113	0.9113	0.9113	0.9144	0.9128	0.9113	0.9707	0.9707	0.9699	0.9715	0.9128	0.9683
0.8906	0.9049	0.9097	0.912	0.9097	0.912	0.912	0.912	0.912	0.9113	0.9105	0.9097	0.9097	0.9675	0.912	0.9144	0.9113
0.8891	0.9065	0.9089	0.9097	0.9105	0.912	0.912	0.9128	0.9128	0.9113	0.9105	0.9097	0.9105	0.9089	0.9089	0.9699	0.9105
0.8938	0.9033	0.9073	0.9073	0.9105	0.9113	0.9136	0.9136	0.9128	0.9128	0.9105	0.9105	0.9105	0.9113	0.9105	0.912	0.9113
0.8986	0.9049	0.9113	0.9128	0.9152	0.9144	0.9128	0.9136	0.916	0.9152	0.9152	0.916	0.9136	0.9128	0.9144	0.9144	0.9136
0.9002	0.912	0.9184	0.9184	0.9208	0.9208	0.9223	0.9223	0.9216	0.9239	0.9231	0.9239	0.9208	0.92	0.9231	0.9216	0.9192
0.9073	0.9255	0.9311	0.9319	0.9319	0.9319	0.9334	0.9326	0.9319	0.9319	0.9326	0.9319	0.9857	0.9319	0.9311	0.9303	0.9303
0.92	0.9834	0.9849	0.9842	0.9857	0.9857	0.9865	0.9865	0.9865	0.9857	0.9857	0.9849	0.9849	0.9849	0.9849	0.9857	0.9857
];


xmin=min(X_Yale);xmax=max(X_Yale);
ymin=min(Y);ymax=max(Y);

[lx,ly,lz] = griddata(X_Yale, Y, Z_Yale, linspace(xmin,xmax)',linspace(ymin,ymax),'v4');


subplot(1, 3, 2)
surf(lx,ly,lz)
title('Acc(\it{s}, \lambda) on ExtYale B')
xlabel('\lambda')
xticks(1:8)
ylabel('\it{s}')
yticks(0:0.2:1)
zlabel('Accuracy')

% UCSD
X_UCSD = [0.1 0.5:0.5:15];
Z_UCSD = [
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



xmin=min(X_UCSD);xmax=max(X_UCSD);
ymin=min(Y);ymax=max(Y);

[lx,ly,lz] = griddata(X_UCSD, Y, Z_UCSD, linspace(xmin,xmax)',linspace(ymin,ymax),'cubic');


subplot(1, 3, 3)
surf(lx,ly,lz)
title('Acc(\it{s}, \lambda) on UCSD')
xlabel('\lambda')
xticks(1:2:15)
ylabel('\it{s}')
yticks(0:0.2:1)
zlabel('Accuracy')

