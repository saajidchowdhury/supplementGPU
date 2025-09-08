clear all;
setenv('TZ', 'America/New_York');
fclose('all');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',50); %get(groot,'factory')
set(groot,'defaultAxesLineWidth',3);
set(groot,'defaultLineLineWidth',3);
set(groot,'defaultLineMarkerSize',50);
set(groot,'defaultErrorbarLineWidth',3);
set(groot,'defaultErrorbarMarkerSize',50);
set(groot,'defaultErrorbarCapSize',20);
set(groot,'defaultAxesView',[0,90]);
set(groot,'defaultAxesBox','on');
set(groot,'defaultTextFontSize',50);
set(groot,'defaultConstantlineLineWidth',3);
set(groot,'defaultConstantlineAlpha',1);
set(groot,'defaultAxesLabelFontSizeMultiplier',1);
set(groot,'defaultFigurePosition',[1282,520,775,686]);
mainhere = string(datetime('now','Format','M-d-y@HH.mm.ss'))+"fig5";
mkdir(mainhere);

load("fig5.mat");

figure; hold on;
xlabel("Number of cores, $n_\textrm{cores}$");
ylabel("Total runtime, $t$ (hours)");
xscale('log');
yscale('log');
xlim([0.9,30]);
ylim([0.02,14]);
xticks([1,10]);
yticks([0.1,1,10]);
plot(ncoress,tcpu/3600);
plot(ncoress,tgpu/3600);
legend("\verb+ode45+","\verb+ode45gpu+");
print(gcf,'-vector','-dsvg',mainhere+"/fig5.svg");
hold off;