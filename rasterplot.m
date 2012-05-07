function rasterplot(time,ch,events)
% RASTERPLOT  plots an image showing time series of binary events
%
% rasterplot(T,CH,EVENTS)
%
% EVENTS is a matrix MxN, for M channels and N samples in time
%









%function raster(time,ch,events)
nTicks = 8;

[I,J] = find(events);
plot(time(J),I,'ko',...
            'MarkerFaceColor','k',...
            'MarkerSize',1.5)

% set ch axis numbering
yInt = floor(length(ch)/nTicks);
yIdx = [(1:yInt:length(ch)-yInt/2) length(ch)];
set(gca,'YTick',yIdx)
set(gca,'YTickLabel',round(ch(yIdx)*10)/10)

axis([time(1) time(end) 1 size(events,1)])


%%%  This method doesn't produce good results for large data matrices - spikes are lost in screen resolution 
% % set extra cells due to refractoriness
% xdel = 1;%max(round(fs*5e-5),1);
% 
% T = (1:size(events,2)).*(1e3/fs);
% CH = (1:size(events,1));
% 
% % find all events
% [ch,t] = find(events);
% for i=1:length(t)
%     events(ch(i), t(i)+1:t(i)+xdel) = 1;
% end
% 
% imagesc(T,CH,events);
% %truesize        % maps display pixel for pixel
% xlabel('Time (ms)')
% ylabel('AN Channel')%'Channel Center Frequency (kHz)')
% 
% colormap([1 1 1; 0 0 0]);       % force black and white
% set(gca,'ydir','normal')        % correct yaxis direction
