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
set(groot,'defaultFigurePosition',[1311,537,746,669]);
mainhere = string(datetime('now','Format','M-d-y@HH.mm.ss'))+"fig9";
mkdir(mainhere);

figure; hold on;
xlabel("Number of trajectories, $n$");
ylabel("Total runtime, $t$ (hours)");
xscale('log');
yscale('log');

load("fig9k80.mat");
tk80 = tendGPU';

load("fig9general.mat");
tgeneral = tendGPU';

plot(ns, tk80(:)/3600);
plot(ns, tgeneral(:)/3600);
xticks([1,100,10000,1000000]);
yticks([0.01,0.1,1,2]);
xlim([0.75,1.5e6]);
ylim([0.0045,2.75]); % ylim([0.0045,2.75]);

legend(["fixed","uniformly sampled"],"Location","Northwest");
print(gcf,'-vector','-dsvg',mainhere+"/fig9.svg");
hold off;