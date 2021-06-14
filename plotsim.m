close all;

datafile = "output.txt";

simfile = importdata(datafile,',',1);
simdata = simfile.data;

%Loop over columns of data
figure;
for n = 2:(size(simdata,2))
    subplot(size(simdata,2),1,n);   %Add a set of axes
    plot(simdata(:,1),simdata(:,n)); %Add the data
    
    ylabel(simfile.colheaders(n));  %Label the Y axis with the column header
    set(gca, 'XScale', 'log')
    xlabel('Frequency (Hz)'); %Label the x axis of the bottom plot
end
