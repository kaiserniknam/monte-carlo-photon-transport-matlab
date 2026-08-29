function [] = Photon_36
% Repository group: 04_compression_beam_and_dpf
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
%%
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: plot mua & mus vs lambda --from selected references for HB and HbO2
% https://omlc.org/spectra/hemoglobin/summary.html
% https://omlc.org/spectra/hemoglobin/
% W. B. Gratzer, Med. Res. Council Labs, Holly Hill, London
% N. Kollias, Wellman Laboratories, Harvard Medical School, Boston
%%
% Martina Meinke, Gerhard Müller, Jürgen Helfmann, and Moritz Friebel, "Empirical model functions to calculate hematocrit-dependent optical properties of human blood," Appl. Opt. 46, 1742-1753 (2007)
% Moritz Friebel and Martina Meinke, "Model function to calculate the refractive index of native hemoglobin in the wavelength range of 250-1100 nm dependent on concentration," Appl. Opt. 45, 2838-2842 (2006)
%%

clc
close all
format long



figure(1)
set_of_eps = calc_epsilon_Hb_HbO2;
eps_Moaveni = epsilon_Hb_HbO2_Moaveni;
eps_Takatani = epsilon_Hb_HbO2_Takatani;
subplot(1,2,1)
semilogy(set_of_eps(:,1),set_of_eps(:,2),'r','DisplayName','HbO2','LineWidth',1.5), hold on
semilogy(set_of_eps(:,1),set_of_eps(:,3),'b','DisplayName','Hb','LineWidth',1.5), hold on
% data - Moaveni
semilogy(eps_Moaveni(:,1),eps_Moaveni(:,2),'sr','DisplayName','HbO2/Moaveni','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k'), hold on
semilogy(eps_Moaveni(:,1),eps_Moaveni(:,3),'sb','DisplayName','Hb/Moaveni','MarkerSize',8,'MarkerFaceColor','b','MarkerEdgeColor','k'), hold on
% data - Takatani
semilogy(eps_Takatani(:,1),eps_Takatani(:,2),'or','DisplayName','HbO2/Takatani','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k'), hold on
semilogy(eps_Takatani(:,1),eps_Takatani(:,3),'ob','DisplayName','Hb/Takatani','MarkerSize',8,'MarkerFaceColor','b','MarkerEdgeColor','k'), hold on
xlabel('wavelength (\lambda, nm)'), ylabel('molar extinction coefficient (\epsilon, cm^{-1}M^{-1})'), title('molar extinction coefficient vs. wavelength for Hb/HbO2')
axis tight, axis square, set(gca,'fontsize',18), axis([200 1000 100 1000000])
legend('show','Location','northeast')
subplot(1,2,2)
semilogy(set_of_eps(:,1),set_of_eps(:,2).*2.303.*150/64500,'r','DisplayName','HbO2','LineWidth',1.5), hold on
semilogy(set_of_eps(:,1),set_of_eps(:,3).*2.303.*150/64500,'b','DisplayName','Hb','LineWidth',1.5), hold on
xlabel('wavelength (\lambda, nm)'), ylabel('absorption coefficient (\mu_a, cm^{-1})'), title('absorption coefficient vs. wavelength for Hb/HbO2')
axis tight, axis square, set(gca,'fontsize',18), axis([200 1000 0 3000])
legend('show','Location','northeast')
clearvars eps_Moaveni eps_Takatani



figure(2)
lmbd = 250:1:1100;
set_of_HCT = [0.84,4.0,5.9,8.6,17.1,25.6,33.2,42.1];
for i_HCT = 1:length(set_of_HCT)
    [mua,mus,g,musp] = calc_muas_based_HCT(set_of_HCT(i_HCT),lmbd);

    subplot(2,2,1)
    semilogy(lmbd,mua/10,'DisplayName',['HCT ',num2str(set_of_HCT(i_HCT)),'%'],'LineWidth',1.5), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('\mu_a, mm^{-1}')
    axis tight, set(gca,'fontsize',14), axis([min(lmbd) max(lmbd) 0.001 1000]),set(gca,'ytick',logspace(-3,+3,7)), set(gca,'xtick',[250 450 650 850 1050]), grid on
    legend('show','Location','northeast','NumColumns',2)

    subplot(2,2,2)
    plot(lmbd,mus/10,'DisplayName',['HCT ',num2str(set_of_HCT(i_HCT)),'%'],'LineWidth',1.5), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('\mu_s, mm^{-1}')
    axis tight, set(gca,'fontsize',14), axis([min(lmbd) max(lmbd) 0 100]),set(gca,'ytick',0:20:100), set(gca,'xtick',[250 450 650 850 1050]), grid on
    legend('show','Location','northeast','NumColumns',2)

    subplot(2,2,3)
    plot(lmbd,g,'DisplayName',['HCT ',num2str(set_of_HCT(i_HCT)),'%'],'LineWidth',1.5), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('g, a.u.')
    axis tight, set(gca,'fontsize',14), axis([min(lmbd) max(lmbd) 0.65 1]),set(gca,'ytick',0.65:0.05:1), set(gca,'xtick',[250 450 650 850 1050]), grid on
    legend('show','Location','southeast','NumColumns',2)

    subplot(2,2,4)
    plot(lmbd,musp/10,'DisplayName',['HCT ',num2str(set_of_HCT(i_HCT)),'%'],'LineWidth',1.5), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('\mu_s'', mm^{-1}')
    axis tight, set(gca,'fontsize',14), axis([min(lmbd) max(lmbd) 0 16]),set(gca,'ytick',0:2:16), set(gca,'xtick',[250 450 650 850 1050]), grid on
    legend('show','Location','northeast','NumColumns',2)

    clearvars mua mus g musp
end
clearvars i_HCT lmbd set_of_HCT



% calculate all optical properties
figure(3)
lmbd = 250:.25:1100;
set_of_HCT = [0.84,4.0,5.9,8.6,17.1,25.6,33.2,42.1];
for i_HCT = 1:length(set_of_HCT)
    [mua,mus,g,musp] = calc_muas_based_HCT(set_of_HCT(i_HCT),lmbd);

    subplot(2,2,1)
    semilogy(lmbd,mua/10,'DisplayName',['HCT ',num2str(set_of_HCT(i_HCT)),'%'],'LineWidth',1.5), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('\mu_a, mm^{-1}')
    axis tight, set(gca,'fontsize',14), axis([min(lmbd) max(lmbd) 0.001 1000]),set(gca,'ytick',logspace(-3,+3,7)), set(gca,'xtick',[250 450 650 850 1050]), grid on
    legend('show','Location','northeast','NumColumns',2)

    subplot(2,2,2)
    plot(lmbd,mus/10,'DisplayName',['HCT ',num2str(set_of_HCT(i_HCT)),'%'],'LineWidth',1.5), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('\mu_s, mm^{-1}')
    axis tight, set(gca,'fontsize',14), axis([min(lmbd) max(lmbd) 0 100]),set(gca,'ytick',0:20:100), set(gca,'xtick',[250 450 650 850 1050]), grid on
    legend('show','Location','northeast','NumColumns',2)

    subplot(2,2,3)
    plot(lmbd,g,'DisplayName',['HCT ',num2str(set_of_HCT(i_HCT)),'%'],'LineWidth',1.5), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('g, a.u.')
    axis tight, set(gca,'fontsize',14), axis([min(lmbd) max(lmbd) 0.65 1]),set(gca,'ytick',0.65:0.05:1), set(gca,'xtick',[250 450 650 850 1050]), grid on
    legend('show','Location','southeast','NumColumns',2)

    % subplot(2,2,4)
    % kapa = lmbd.*(1e-9).*mua.*(1e+2)/4/pi;
    % rfindx = nan(size(lmbd));
    % for i_lmbda = 1:length(lmbd)
    %     s = (2/pi).*(kapa./lmbd).*(1./(1-(lmbd./lmbd(i_lmbda)).^2)).*mean(diff(lmbd));
    %     s = s(~isinf(s));
    %     rfindx(i_lmbda) = 1+sum(s); clearvars s
    % end
    % plot(lmbd,rfindx,'DisplayName',['HCT ',num2str(set_of_HCT(i_HCT)),'%'],'LineWidth',1.5), hold on
    % xlabel('wavelength (\lambda, nm)'), ylabel('n from kk, a.u.')
    % axis tight, set(gca,'fontsize',14), % axis([min(lmbd) max(lmbd) 1.30 1.50]),set(gca,'ytick',1.30:0.02:1.50), set(gca,'xtick',[250 450 650 850 1050]), grid on
    % legend('show','Location','southeast','NumColumns',2)

    subplot(2,2,4)
    plot(lmbd,refractive_index_water(lmbd).*((Specific_Refractive_Increment_beta(lmbd).*set_of_HCT(i_HCT)./3)+1),'DisplayName',['HCT ',num2str(set_of_HCT(i_HCT)),'%'],'LineWidth',1.5), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('n, a.u.')
    axis tight, set(gca,'fontsize',14), axis([min(lmbd) max(lmbd) 1.30 1.50]),set(gca,'ytick',1.30:0.02:1.50), set(gca,'xtick',[250 450 650 850 1050]), grid on
    legend('show','Location','northeast','NumColumns',2)

    clearvars mua mus g musp rfindx kapa i_lmbda
end
clearvars i_HCT lmbd set_of_HCT



% calculate n (refractive index)
figure(4)
lmbd = 250:.25:1100;
set_of_Hb = [0, 4.6, 10.4, 16.5, 28.7];
for i_Hb = 1:length(set_of_Hb)
    plot(lmbd,refractive_index_water_two(lmbd).*(Specific_Refractive_Increment_beta(lmbd).*set_of_Hb(i_Hb)+1),'DisplayName',['Hb ',num2str(set_of_Hb(i_Hb)),' g/dL'],'LineWidth',1.5), hold on
    xlabel('wavelength (\lambda, nm)'), ylabel('n, a.u.')
    axis tight, set(gca,'fontsize',14), axis([min(lmbd) max(lmbd) 1.30 1.50]),set(gca,'ytick',1.30:0.02:1.50), set(gca,'xtick',[250 450 650 850 1050]), grid on
    legend('show','Location','southeast','NumColumns',2)

    clearvars mua mus g musp rfindx kapa i_lmbda
end
clearvars i_Hb lmbd set_of_Hb


% comparison
figure(5)
lmbd = 250:.25:1100;
set_of_HCT = [42.1];
set_of_eps = calc_epsilon_Hb_HbO2;
semilogy(set_of_eps(:,1),set_of_eps(:,2).*2.303.*150/64500,'r','DisplayName','Data: HbO2','LineWidth',1.5), hold on
semilogy(set_of_eps(:,1),set_of_eps(:,3).*2.303.*150/64500,'b','DisplayName','Data: Hb','LineWidth',1.5), hold on
for i_HCT = 1:length(set_of_HCT)
    [mua,~,~,~] = calc_muas_based_HCT(set_of_HCT(i_HCT),lmbd);
    semilogy(lmbd,mua,'DisplayName',['HCT ',num2str(set_of_HCT(i_HCT)),'%'],'LineWidth',1.5,'Color','g'), hold on
    clearvars mua mus g musp rfindx kapa i_lmbda
end
clearvars i_HCT lmbd set_of_HCT
xlabel('wavelength (\lambda, nm)'), ylabel('absorption coefficient (\mu_a, cm^{-1})'), title('absorption coefficient vs. wavelength for Hb/HbO2')
axis tight, axis square, set(gca,'fontsize',18), axis([200 1000 0 3000]), grid on
legend('show','Location','northeast')
clearvars eps_Moaveni eps_Takatani
end

function [out] = calc_epsilon_Hb_HbO2()
out = [...
250	106112	112736
252	105552	112736
254	107660	112736
256	109788	113824
258	112944	115040
260	116376	116296
262	120188	117564
264	124412	118876
266	128696	120208
268	133064	121544
270	136068	122880
272	137232	123096
274	138408	121952
276	137424	120808
278	135820	119840
280	131936	118872
282	127720	117628
284	122280	114820
286	116508	112008
288	108484	107140
290	104752	98364
292	98936	91636
294	88136	85820
296	79316	77100
298	70884	69444
300	65972	64440
302	63208	61300
304	61952	58828
306	62352	56908
308	62856	57620
310	63352	59156
312	65972	62248
314	69016	65344
316	72404	68312
318	75536	71208
320	78752	74508
322	82256	78284
324	85972	82060
326	89796	85592
328	93768	88516
330	97512	90856
332	100964	93192
334	103504	95532
336	104968	99792
338	106452	104476
340	107884	108472
342	109060	110996
344	110092	113524
346	109032	116052
348	107984	118752
350	106576	122092
352	105040	125436
354	103696	128776
356	101568	132120
358	97828	133632
360	94744	134940
362	92248	136044
364	89836	136972
366	88484	137900
368	87512	138856
370	88176	139968
372	91592	141084
374	95140	142196
376	98936	143312
378	103432	144424
380	109564	145232
382	116968	145232
384	125420	148668
386	135132	153908
388	148100	159544
390	167748	167780
392	189740	180004
394	212060	191540
396	231612	202124
398	248404	212712
400	266232	223296
402	284224	236188
404	308716	253368
406	354208	270548
408	422320	287356
410	466840	303956
412	500200	321344
414	524280	342596
416	521880	363848
418	515520	385680
420	480360	407560
422	431880	429880
424	376236	461200
426	326032	481840
428	283112	500840
430	246072	528600
432	214120	552160
434	165332	552160
436	132820	547040
438	119140	501560
440	102580	413280
442	92780	363240
444	81444	282724
446	76324	237224
448	67044	173320
450	62816	103292
452	58864	62640
454	53552	36170
456	49496	30698.8
458	47496	25886.4
460	44480	23388.8
462	41320	20891.2
464	39807.2	19260.8
466	37073.2	18142.4
468	34870.8	17025.6
470	33209.2	16156.4
472	31620	15310
474	30113.6	15048.4
476	28850.8	14792.8
478	27718	14657.2
480	26629.2	14550
482	25701.6	14881.2
484	25180.4	15212.4
486	24669.6	15543.6
488	24174.8	15898
490	23684.4	16684
492	23086.8	17469.6
494	22457.6	18255.6
496	21850.4	19041.2
498	21260	19891.2
500	20932.8	20862
502	20596.4	21832.8
504	20418	22803.6
506	19946	23774.4
508	19996	24745.2
510	20035.2	25773.6
512	20150.4	26936.8
514	20429.2	28100
516	21001.6	29263.2
518	22509.6	30426.4
520	24202.4	31589.6
522	26450.4	32851.2
524	29269.2	34397.6
526	32496.4	35944
528	35990	37490
530	39956.8	39036.4
532	43876	40584
534	46924	42088
536	49752	43592
538	51712	45092
540	53236	46592
542	53292	48148
544	52096	49708
546	49868	51268
548	46660	52496
550	43016	53412
552	39675.2	54080
554	36815.2	54520
556	34476.8	54540
558	33456	54164
560	32613.2	53788
562	32620	52276
564	33915.6	50572
566	36495.2	48828
568	40172	46948
570	44496	45072
572	49172	43340
574	53308	41716
576	55540	40092
578	54728	38467.6
580	50104	37020
582	43304	35676.4
584	34639.6	34332.8
586	26600.4	32851.6
588	19763.2	31075.2
590	14400.8	28324.4
592	10468.4	25470
594	7678.8	22574.8
596	5683.6	19800
598	4504.4	17058.4
600	3200	14677.2
602	2664	13622.4
604	2128	12567.6
606	1789.2	11513.2
608	1647.6	10477.6
610	1506	9443.6
612	1364.4	8591.2
614	1222.8	7762
616	1110	7344.8
618	1026	6927.2
620	942	6509.6
622	858	6193.2
624	774	5906.8
626	707.6	5620
628	658.8	5366.8
630	610	5148.8
632	561.2	4930.8
634	512.4	4730.8
636	478.8	4602.4
638	460.4	4473.6
640	442	4345.2
642	423.6	4216.8
644	405.2	4088.4
646	390.4	3965.08
648	379.2	3857.6
650	368	3750.12
652	356.8	3642.64
654	345.6	3535.16
656	335.2	3427.68
658	325.6	3320.2
660	319.6	3226.56
662	314	3140.28
664	308.4	3053.96
666	302.8	2967.68
668	298	2881.4
670	294	2795.12
672	290	2708.84
674	285.6	2627.64
676	282	2554.4
678	279.2	2481.16
680	277.6	2407.92
682	276	2334.68
684	274.4	2261.48
686	272.8	2188.24
688	274.4	2115
690	276	2051.96
692	277.6	2000.48
694	279.2	1949.04
696	282	1897.56
698	286	1846.08
700	290	1794.28
702	294	1741
704	298	1687.76
706	302.8	1634.48
708	308.4	1583.52
710	314	1540.48
712	319.6	1497.4
714	325.2	1454.36
716	332	1411.32
718	340	1368.28
720	348	1325.88
722	356	1285.16
724	364	1244.44
726	372.4	1203.68
728	381.2	1152.8
730	390	1102.2
732	398.8	1102.2
734	407.6	1102.2
736	418.8	1101.76
738	432.4	1100.48
740	446	1115.88
742	459.6	1161.64
744	473.2	1207.4
746	487.6	1266.04
748	502.8	1333.24
750	518	1405.24
752	533.2	1515.32
754	548.4	1541.76
756	562	1560.48
758	574	1560.48
760	586	1548.52
762	598	1508.44
764	610	1459.56
766	622.8	1410.52
768	636.4	1361.32
770	650	1311.88
772	663.6	1262.44
774	677.2	1213
776	689.2	1163.56
778	699.6	1114.8
780	710	1075.44
782	720.4	1036.08
784	730.8	996.72
786	740	957.36
788	748	921.8
790	756	890.8
792	764	859.8
794	772	828.8
796	786.4	802.96
798	807.2	782.36
800	816	761.72
802	828	743.84
804	836	737.08
806	844	730.28
808	856	723.52
810	864	717.08
812	872	711.84
814	880	706.6
816	887.2	701.32
818	901.6	696.08
820	916	693.76
822	930.4	693.6
824	944.8	693.48
826	956.4	693.32
828	965.2	693.2
830	974	693.04
832	982.8	692.92
834	991.6	692.76
836	1001.2	692.64
838	1011.6	692.48
840	1022	692.36
842	1032.4	692.2
844	1042.8	691.96
846	1050	691.76
848	1054	691.52
850	1058	691.32
852	1062	691.08
854	1066	690.88
856	1072.8	690.64
858	1082.4	692.44
860	1092	694.32
862	1101.6	696.2
864	1111.2	698.04
866	1118.4	699.92
868	1123.2	701.8
870	1128	705.84
872	1132.8	709.96
874	1137.6	714.08
876	1142.8	718.2
878	1148.4	722.32
880	1154	726.44
882	1159.6	729.84
884	1165.2	733.2
886	1170	736.6
888	1174	739.96
890	1178	743.6
892	1182	747.24
894	1186	750.88
896	1190	754.52
898	1194	758.16
900	1198	761.84
902	1202	765.04
904	1206	767.44
906	1209.2	769.8
908	1211.6	772.16
910	1214	774.56
912	1216.4	776.92
914	1218.8	778.4
916	1220.8	778.04
918	1222.4	777.72
920	1224	777.36
922	1225.6	777.04
924	1227.2	776.64
926	1226.8	772.36
928	1224.4	768.08
930	1222	763.84
932	1219.6	752.28
934	1217.2	737.56
936	1215.6	722.88
938	1214.8	708.16
940	1214	693.44
942	1213.2	678.72
944	1212.4	660.52
946	1210.4	641.08
948	1207.2	621.64
950	1204	602.24
952	1200.8	583.4
954	1197.6	568.92
956	1194	554.48
958	1190	540.04
960	1186	525.56
962	1182	511.12
964	1178	495.36
966	1173.2	473.32
968	1167.6	451.32
970	1162	429.32
972	1156.4	415.28
974	1150.8	402.28
976	1144	389.288
978	1136	374.944
980	1128	359.656
982	1120	344.372
984	1112	329.084
986	1102.4	313.796
988	1091.2	298.508
990	1080	283.22
992	1068.8	267.932
994	1057.6	252.648
996	1046.4	237.36
998	1035.2	222.072
1000	1024	206.784];
end
function [out] = epsilon_Hb_HbO2_Moaveni()
% J.M. Schmitt, "Optical Measurement of Blood Oxygenation by Implantable Telemetry," Technical Report G558-15, Stanford."
% M.K. Moaveni, "A Multiple Scattering Field Theory Applied to Whole Blood," Ph.D. dissertation, Dept. of Electrical Engineering, University of Washington, 1970.
out = [...
630	680	4280
640	440	3640
650	380	3420
660	320	3200
670	320	3080
680	320	2960
690	280	2560
700	320	2160
710	340	1840
720	360	1520
730	400	1500
740	440	1520
750	520	1620
760	600	1720
770	660	1420
780	720	1120
790	760	1020
800	800	920
810	860	880
820	920	840
830	980	840
840	1040	840
850	1060	800
860	1080	840
870	1120	840
880	1160	840
890	1180	860
900	1200	880
910	1220	920
920	1240	880
930	1240	800
940	1200	800
950	1200	720];
end
function [out] = epsilon_Hb_HbO2_Takatani()
% S. Takatani and M. D. Graham, "Theoretical analysis of diffuse reflectance from a two-layer tissue model," IEEE Trans. Biomed. Eng., BME-26, 656--664, (1987).
out = [...
450	68000	58000
460	45040	20600
480	27360	13360
500	20200	16360
507	19240	19240
510	19040	20000
520	23520	25080
522	25680	25680
540	57080	41120
542	57480	44000
549	49840	49840
555	36000	52160
560	33880	50160
569	45080	45080
577	61480	36800
579	54920	35440
586	28920	28920
600	3200	14600
605	1860	9496
615	1152	5776
625	732	4400
635	488	3796
645	396	3436
655	340	3244
665	292	3156
675	288	3028
685	272	2796
695	280	2424
705	300	1988
715	328	1628
725	368	1464
735	412	1464
745	480	1616
755	556	1756
765	616	1640
775	684	1340
785	736	1040
795	776	964
805	880	896
815	880	880
825	952	832
835	996	820
845	1048	820
855	1068	820
865	1116	820
875	1140	848
885	1168	832
895	1188	884
905	1208	896
915	1220	924
925	1228	860
935	1216	848
945	1212	756
955	1196	704
965	1176	616
975	1148	552
985	1108	424
995	1052	372];
end
function [mua_st,mus_st,g_st,musp_st] = Data_of_Standard_Optical_Parameters (lmbda)
% Martina Meinke, Gerhard Müller, Jürgen Helfmann, and Moritz Friebel, "Empirical model functions to calculate hematocrit-dependent optical properties of human blood," Appl. Opt. 46, 1742-1753 (2007)
db = [
250	7.52	27.9	0.877	3.43
255	7.93	27.2	0.88	3.27
260	8.6	26.6	0.883	3.1
265	9.41	25.7	0.883	3.01
270	9.82	25.4	0.882	2.99
275	9.93	25.2	0.884	2.93
280	9.74	25.3	0.889	2.8
285	9.19	25.5	0.899	2.59
290	8.24	25.9	0.912	2.28
295	7.04	26.7	0.9268	1.96
300	6	27.4	0.9389	1.68
305	5.52	27.7	0.9449	1.53
310	5.56	27.5	0.9442	1.54
315	5.95	26.9	0.9398	1.62
320	6.5	26	0.9335	1.73
325	7.13	25.4	0.928	1.83
330	7.69	25.2	0.9239	1.92
335	8.1	24.9	0.9208	1.97
340	8.38	24.7	0.9192	1.99
345	8.48	24.6	0.9191	1.99
350	8.37	24.6	0.9215	1.93
355	8.07	24.7	0.9257	1.84
360	7.69	24.7	0.9308	1.71
365	7.4	24.7	0.9345	1.61
370	7.37	24.5	0.9356	1.58
375	7.72	24	0.9323	1.62
380	8.58	23.2	0.9242	1.76
385	9.98	22.2	0.9105	1.99
390	11.86	21.1	0.892	2.29
395	14.21	20.1	0.869	2.64
400	16.86	19.1	0.845	2.96
405	19.66	18.3	0.824	3.23
410	21.8	17.9	0.81	3.39
415	22.66	18.6	0.812	3.48
420	21.54	19.8	0.831	3.35
425	18.62	21.8	0.861	3.03
430	14.93	23.8	0.893	2.55
435	11.52	25.3	0.92	2.02
440	8.95	26.8	0.9405	1.6
445	7.11	28	0.9543	1.28
450	5.79	29	0.9629	1.08
455	4.83	30	0.9685	0.946
460	4.11	30.5	0.9724	0.841
465	3.56	31	0.9756	0.757
470	3.13	31.3	0.9777	0.699
475	2.79	31.5	0.9793	0.652
480	2.52	31.8	0.9806	0.619
490	2.2	32	0.982	0.575
500	2.04	32.2	0.9831	0.545
510	1.98	32.2	0.9835	0.531
520	2.51	31.5	0.9838	0.51
530	3.97	29.5	0.9794	0.609
540	5.05	28.2	0.9755	0.691
550	4.4	29	0.9779	0.642
560	3.6	30.1	0.9804	0.59
570	4.56	28.8	0.9777	0.641
580	4.56	28.9	0.9771	0.662
590	1.9	32.1	0.9827	0.556
600	0.478	33.9	0.9854	0.496
610	0.17	34.2	0.9858	0.487
620	0.0812	34.3	0.9861	0.477
630	0.0496	34.2	0.9863	0.469
640	0.0348	34.3	0.9865	0.462
660	0.0251	34.1	0.9868	0.45
670	0.0239	33.9	0.9871	0.439
690	0.0243	33.6	0.9872	0.43
700	0.0246	33.4	0.9872	0.427
720	0.0284	32.9	0.9871	0.426
740	0.036	32.4	0.987	0.422
760	0.0461	31.8	0.9868	0.42
780	0.0558	31.2	0.9867	0.415
800	0.0641	30.8	0.9868	0.407
820	0.0762	30.5	0.9867	0.407
840	0.0853	30.5	0.9867	0.405
860	0.0953	29.7	0.9864	0.403
880	0.104	29.8	0.9861	0.413
900	0.106	28.4	0.9857	0.405
920	0.111	28.1	0.9854	0.411
940	0.117	27.7	0.9847	0.423
960	0.125	26.8	0.984	0.429
980	0.133	26.1	0.9836	0.428
1000	0.128	25.8	0.9837	0.421
1020	0.118	25.7	0.9839	0.412
1040	0.104	25.3	0.9837	0.412
1060	0.0879	25	0.9838	0.406
1080	0.0795	24.6	0.9841	0.392
1100	0.074	24.5	0.9842	0.387
];
mua_st  = interp1(db(:,1),db(:,2),lmbda,'linear','extrap'); mua_st  =  mua_st*10;
mus_st  = interp1(db(:,1),db(:,3),lmbda,'linear','extrap'); mus_st  =  mus_st*10;
g_st    = interp1(db(:,1),db(:,4),lmbda,'linear','extrap'); g_st    =       g_st;
musp_st = interp1(db(:,1),db(:,5),lmbda,'linear','extrap'); musp_st = musp_st*10;
end
function [mua,mus,g,musp] = calc_muas_based_HCT(HCT,lambda)
% Martina Meinke, Gerhard Müller, Jürgen Helfmann, and Moritz Friebel, "Empirical model functions to calculate hematocrit-dependent optical properties of human blood," Appl. Opt. 46, 1742-1753 (2007)
[mua_st,mus_st,g_st,musp_st] = Data_of_Standard_Optical_Parameters (lambda);
mua = nan(size(lambda));
mus = nan(size(lambda));
g   = nan(size(lambda));
musp= nan(size(lambda));
for i_lambda = 1:length(lambda)
    % mu_a
    if     ((250<=lambda(i_lambda)&&lambda(i_lambda)<=400)||(430<=lambda(i_lambda)&&lambda(i_lambda)<=600))&&(0.84<=HCT&&HCT<=42.1)
        mua(i_lambda) = 0.1233.*mua_st(i_lambda).*HCT;
    elseif  (600<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(0.84<=HCT&&HCT<=42.1)
        mua(i_lambda) = 0.1206.*mua_st(i_lambda).*HCT;
    elseif  (400< lambda(i_lambda)&&lambda(i_lambda)< 430 )&&(0.84<=HCT&&HCT<=42.1)
        mua(i_lambda) = mean([0.1233 0.1206]).*mua_st(i_lambda).*HCT;
    else
        mua(i_lambda) = nan;
    end
    % mu_s
    if      (250<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(0.84<=HCT&&HCT<=17.1)
        mus(i_lambda) = (-0.0015.*HCT.^2+0.1268.*HCT).*mus_st(i_lambda);
    elseif  (600<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(17.1<=HCT&&HCT<=42.1)
        mus(i_lambda) = (-0.0015.*HCT.^2+0.1268.*HCT).*mus_st(i_lambda);
    elseif  (250<=lambda(i_lambda)&&lambda(i_lambda)< 600 )&&(17.1<=HCT&&HCT<=42.1)
        mus(i_lambda) = (-0.0015.*HCT.^2+0.1268.*HCT).*mus_st(i_lambda);
    else
        mus(i_lambda) = nan;
    end
    % mu_sp
    if      (250<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(0.84<=HCT&&HCT<=17.1)
        musp(i_lambda) = 0.1167.*HCT.*musp_st(i_lambda);
    elseif  (600<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(17.1<=HCT&&HCT<=42.1)
        musp(i_lambda) = 0.1167.*HCT.*musp_st(i_lambda);
    elseif  (250<=lambda(i_lambda)&&lambda(i_lambda)< 600 )&&(17.1<=HCT&&HCT<=42.1)
        musp(i_lambda) = 0.1167.*HCT.*musp_st(i_lambda);
    else
        musp(i_lambda) = nan;
    end
    % g
    g = (1-musp./mus);
    % if      (600<=lambda(i_lambda)&&lambda(i_lambda)<=1100)&&(0.84<=HCT&&HCT<=42.1)
    %     g(i_lambda) = (((-2.684e-6).*HCT.^2)+((-2.373e-4).*HCT)+1.003).*g_st(i_lambda);
    % elseif  (250<=lambda(i_lambda)&&lambda(i_lambda)< 600 )&&(0.84<=HCT&&HCT<=42.1)
    %     g(i_lambda) = (((-2.684e-6).*HCT.^2)+((-2.373e-4).*HCT)+1.003).*g_st(i_lambda);
    % else
    %     g(i_lambda) = nan;
    % end
end
end
function [beta_st] = Specific_Refractive_Increment_beta(lmbda)
% Moritz Friebel and Martina Meinke, "Model function to calculate the refractive index of native hemoglobin in the wavelength range of 250-1100 nm dependent on concentration," Appl. Opt. 45, 2838-2842 (2006)
db = [...
250	0.00221
255	0.002155
260	0.002105
265	0.002069
270	0.002048
275	0.002042
280	0.002044
285	0.002047
290	0.002047
295	0.002037
300	0.00202
305	0.001999
310	0.001998
320	0.002007
330	0.002021
340	0.00201
350	0.001989
355	0.001985
360	0.001983
365	0.001912
370	0.00186
375	0.001816
380	0.001774
385	0.001732
390	0.001694
395	0.001668
400	0.001664
405	0.001701
410	0.001799
415	0.001985
420	0.002117
425	0.002195
430	0.002273
435	0.002227
440	0.00221
445	0.002184
450	0.002156
455	0.002131
460	0.002109
465	0.002092
470	0.002078
475	0.002067
480	0.002056
485	0.002045
490	0.002033
495	0.002019
500	0.002005
510	0.002009
520	0.001983
530	0.001966
540	0.001981
550	0.001998
560	0.001992
570	0.001988
580	0.002004
590	0.002015
600	0.001988
610	0.001967
620	0.001964
630	0.00196
640	0.001954
660	0.001958
680	0.00197
700	0.001992
720	0.001979
740	0.001955
760	0.001958
780	0.00196
800	0.001939
820	0.00192
840	0.001935
860	0.001951
880	0.001982
900	0.001998
920	0.002011
940	0.002015
960	0.002021
980	0.002017
1000	0.002052
1020	0.002049
1040	0.002044
1060	0.00204
1080	0.002044
1100	0.002056
];
beta_st = interp1(db(:,1),db(:,2),lmbda,'linear','extrap');
end
function [n_water] = refractive_index_water(lambda_nm)
% Computes the refractive index of water as a function of wavelength in nm
% George M. Hale and Marvin R. Querry, "Optical Constants of Water in the 200-nm to 200-μm Wavelength Region," Appl. Opt. 12, 555-563 (1973)
% from https://refractiveindex.info/?book=H2O&page=Hale&shelf=main&utm_source=chatgpt.com
db = [...
    0.200	1.396
    0.225	1.373
    0.250	1.362
    0.275	1.354
    0.300	1.349
    0.325	1.346
    0.350	1.343
    0.375	1.341
    0.400	1.339
    0.425	1.338
    0.450	1.337
    0.475	1.336
    0.500	1.335
    0.525	1.334
    0.550	1.333
    0.575	1.333
    0.600	1.332
    0.625	1.332
    0.650	1.331
    0.675	1.331
    0.700	1.331
    0.725	1.330
    0.750	1.330
    0.775	1.330
    0.800	1.329
    0.825	1.329
    0.850	1.329
    0.875	1.328
    0.900	1.328
    0.925	1.328
    0.950	1.327
    0.975	1.327
    1.0	1.327
    1.2	1.324
    1.4	1.321
    1.6	1.317
    1.8	1.312
    2.0	1.306
    2.2	1.296
    2.4	1.279
    2.6	1.242
    2.65	1.219
    2.70	1.188
    2.75	1.157
    2.80	1.142
    2.85	1.149
    2.90	1.201
    2.95	1.292
    3.00	1.371
    3.05	1.426
    3.10	1.467
    3.15	1.483
    3.20	1.478
    3.25	1.467
    3.30	1.450
    3.35	1.432
    3.40	1.420
    3.45	1.410
    3.50	1.400
    3.6	1.385
    3.7	1.374
    3.8	1.364
    3.9	1.357
    4.0	1.351
    4.1	1.346
    4.2	1.342
    4.3	1.338
    4.4	1.334
    4.5	1.332
    4.6	1.330
    4.7	1.330
    4.8	1.330
    4.9	1.328
    5.0	1.325
    5.1	1.322
    5.2	1.317
    5.3	1.312
    5.4	1.305
    5.5	1.298
    5.6	1.289
    5.7	1.277
    5.8	1.262
    5.9	1.248
    6.0	1.265
    6.1	1.319
    6.2	1.363
    6.3	1.357
    6.4	1.347
    6.5	1.339
    6.6	1.334
    6.7	1.329
    6.8	1.324
    6.9	1.321
    7.0	1.317
    7.1	1.314
    7.2	1.312
    7.3	1.309
    7.4	1.307
    7.5	1.304
    7.6	1.302
    7.7	1.299
    7.8	1.297
    7.9	1.294
    8.0	1.291
    8.2	1.286
    8.4	1.281
    8.6	1.275
    8.8	1.269
    9.0	1.262
    9.2	1.255
    9.4	1.247
    9.6	1.239
    9.8	1.229
    10.0	1.218
    10.5	1.185
    11.0	1.153
    11.5	1.126
    12.0	1.111
    12.5	1.123
    13.0	1.146
    13.5	1.177
    14.0	1.210
    14.5	1.241
    15.0	1.270
    15.5	1.297
    16.0	1.325
    16.5	1.351
    17.0	1.376
    17.5	1.401
    18.0	1.423
    18.5	1.443
    19.0	1.461
    19.5	1.476
    20.0	1.480
    21.0	1.487
    22	1.500
    23	1.511
    24	1.521
    25	1.531
    26	1.539
    27	1.545
    28	1.549
    29	1.551
    30	1.551
    32	1.546
    34	1.536
    36	1.527
    38	1.522
    40	1.519
    42	1.522
    44	1.530
    46	1.541
    48	1.555
    50	1.587
    60	1.703
    70	1.821
    80	1.886
    90	1.924
    100	1.957
    110	1.966
    120	2.004
    130	2.036
    140	2.056
    150	2.069
    160	2.081
    170	2.094
    180	2.107
    190	2.119
    200	2.130
    ];
n_water = interp1(db(:,1).*1000,db(:,2),lambda_nm,'linear','extrap');
end
function [n_water] = refractive_index_water_two(lambda_nm)
% Computes the refractive index of water as a function of wavelength in nm
% George M. Hale and Marvin R. Querry, "Optical Constants of Water in the 200-nm to 200-μm Wavelength Region," Appl. Opt. 12, 555-563 (1973)
db = [...
    0.2     1.1e-7   1.396
    0.225   4.9e-8   1.373
    0.25    3.35e-8  1.362
    0.275   2.35e-8  1.354
    0.3     1.6e-8   1.349
    0.325   1.08e-8  1.346
    0.35    6.5e-9   1.343
    0.375   3.5e-9   1.341
    0.4     1.86e-9  1.339
    0.425   1.3e-9   1.338
    0.45    1.02e-9  1.337
    0.475   9.35e-10 1.336
    0.5     1.00e-9  1.335
    0.525   1.32e-9  1.334
    0.55    1.96e-9  1.333
    0.575   3.60e-9  1.333
    0.6     1.09e-8  1.332
    0.625   1.39e-8  1.332
    0.65    1.64e-8  1.331
    0.675   2.23e-8  1.331
    0.7     3.35e-8  1.331
    0.725   9.15e-8  1.33
    0.75    1.56e-7  1.33
    0.775   1.48e-7  1.33
    0.8     1.25e-7  1.329
    0.825   1.82e-7  1.329
    0.85    2.93e-7  1.329
    0.875   3.91e-7  1.328
    0.9     4.86e-7  1.328
    0.925   1.06e-6  1.328
    0.95    2.93e-6  1.327
    0.975   3.48e-6  1.327
    1       2.89e-6  1.327
    1.2     9.89e-6  1.324
    1.4     1.38e-4  1.321
    1.6     8.55e-5  1.317
    1.8     1.15e-4  1.312
    2       1.1e-3   1.306
    2.2     2.89e-4  1.296
    2.4     9.56e-4  1.279
    2.6     3.17e-3  1.242
    2.65    6.7e-3   1.219
    2.7     0.019    1.188
    2.75    0.059    1.157
    2.8     0.115    1.142
    2.85    0.185    1.149
    2.9     0.268    1.201
    2.95    0.298    1.292
    3       0.272    1.371
    3.05    0.24     1.426
    3.1     0.192    1.467
    3.15    0.135    1.483
    3.2     0.0924   1.478
    3.25    0.061    1.467
    3.3     0.0368   1.45
    3.35    0.0261   1.432
    3.4     0.0195   1.42
    3.45    0.0132   1.41
    3.5     0.0094   1.4
    3.6     0.00515  1.385
    3.7     0.0036   1.374
    3.8     0.0034   1.364
    3.9     0.0038   1.357
    4       0.0046   1.351
    4.1     0.00562  1.346
    4.2     0.00688  1.342
    4.3     0.00845  1.338
    4.4     0.0103   1.334
    4.5     0.0134   1.332
    4.6     0.0147   1.33
    4.7     0.0157   1.33
    4.8     0.015    1.33
    4.9     0.0137   1.328
    5       0.0124   1.325
    5.1     0.0111   1.322
    5.2     0.0101   1.317
    5.3     0.0098   1.312
    5.4     0.0103   1.305
    5.5     0.0116   1.298
    5.6     0.0142   1.289
    5.7     0.0203   1.277
    5.8     0.033    1.262
    5.9     0.0622   1.248
    6       0.107    1.265
    6.1     0.131    1.319
    6.2     0.088    1.363
    6.3     0.057    1.357
    6.4     0.0449   1.347
    6.5     0.0392   1.339
    6.6     0.0356   1.334
    6.7     0.0337   1.329
    6.8     0.0327   1.324
    6.9     0.0322   1.321
    7       0.032    1.317
    7.1     0.032    1.314
    7.2     0.0321   1.312
    7.3     0.0322   1.309
    7.4     0.0324   1.307
    7.5     0.0326   1.304
    7.6     0.0328   1.302
    7.7     0.0331   1.299
    7.8     0.0335   1.297
    7.9     0.0339   1.294
    8       0.0343   1.291
    8.2     0.0351   1.286
    8.4     0.0361   1.281
    8.6     0.0372   1.275
    8.8     0.0385   1.269
    9       0.0399   1.262
    9.2     0.0415   1.255
    9.4     0.0433   1.247
    9.6     0.0454   1.239
    9.8     0.0479   1.229
    10      0.0508   1.218
    10.5    0.0662   1.185
    11      0.0968   1.153
    11.5    0.142    1.126
    12      0.199    1.111
    12.5    0.259    1.123
    13      0.305    1.146
    13.5    0.343    1.177
    14      0.37     1.21
    14.5    0.388    1.241
    15      0.402    1.27
    15.5    0.414    1.297
    16      0.422    1.325
    16.5    0.428    1.351
    17      0.429    1.376
    17.5    0.429    1.401
    18      0.426    1.423
    18.5    0.421    1.443
    19      0.414    1.461
    19.5    0.404    1.476
    20      0.393    1.48
    21      0.382    1.487
    22      0.373    1.5
    23      0.367    1.511
    24      0.361    1.521
    25      0.356    1.531
    26      0.35     1.539
    27      0.344    1.545
    28      0.338    1.549
    29      0.333    1.551
    30      0.328    1.551
    32      0.324    1.546
    34      0.329    1.536
    36      0.343    1.527
    38      0.361    1.522
    40      0.385    1.519
    42      0.409    1.522
    44      0.436    1.53
    46      0.462    1.541
    48      0.488    1.555
    50      0.514    1.587
    60      0.587    1.703
    70      0.576    1.821
    80      0.547    1.886
    90      0.536    1.924
    100     0.532    1.957
    110     0.531    1.966
    120     0.526    2.004
    130     0.514    2.036
    140     0.5      2.056
    150     0.495    2.069
    160     0.496    2.081
    170     0.497    2.094
    180     0.499    2.107
    190     0.501    2.119
    200     0.504    2.13...
    ];
n_water = interp1(db(:,1).*1000,db(:,3),lambda_nm,'linear','extrap');
end
