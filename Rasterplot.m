function [figH] = Rasterplot(Cell,param,color,fs)

% This function plots a rasterplot of the neurons
%
% Boles?aw Osi?ski (2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%
% Cell      - Group of cell where we're taking the activity;
% param     - set of network parameters
% color     - [r,g,b] numerical color values
% fs        - fontsize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scrsz = get(0,'ScreenSize');
figH = figure;
set(figH,'position',[0,400,scrsz(3)-0.6*scrsz(3),scrsz(4)-0.8*scrsz(4)]);
cla;
hold on;
figtitle = [Cell{1}.label,' cells'];
title(figtitle,'fontsize',16);
for ii = 1:length(Cell)
    J = find(Cell{ii}.S);
    for jj = 1:length(J)
        spkx = [J(jj),J(jj)] .* param.dt;
        spky = [ii,ii + 0.9];
        line(spkx,spky,'color',color,'LineWidth',1);
    end
end

axis([0,param.tsim + param.dt,0,length(Cell) + 2]);
set(gca,'fontsize',fs)
xlabel('time (ms)','fontsize',fs);ylabel('neuron','fontsize',fs);
