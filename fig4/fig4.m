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
set(groot,'defaultFigurePosition',[1214,380,843,826]);
mainhere = string(datetime('now','Format','M-d-y@HH.mm.ss'))+"fig4";
mkdir(mainhere);

figure; hold on;
xlabel("Number of trajectories, $n$");
ylabel("Cumulative runtime, $t$ (hours)");

load("fig4olderversion.mat");
tcpu = tendCPU';
tolderversion = tendGPU';

load("fig4versiongpu.mat");
tversiongpu = tendGPU';

plot(1:ntheta*nphi, cumsum(tcpu(:))/3600);
plot(1:ntheta*nphi, cumsum(tolderversion(:))/3600);
plot(1:ntheta*nphi, cumsum(tversiongpu(:))/3600);
xlim([-50,2530]);
ylim([-0.25,11]);

legend(["\verb+ode45+","older version","\verb+ode45gpu+"],"Location","Northwest");
print(gcf,'-vector','-dsvg',mainhere+"/time.svg");
hold off;