function recording_report(params,bidsID)

import mlreportgen.report.*
import mlreportgen.dom.*

f = strcat('report-' ,bidsID);
rpt = Report(strcat('/tmp/', f), 'pdf');

% Title page
tp = TitlePage;
tp.Title = ['Report of ' bidsID];
append(rpt,tp);

%% Preprocessing report
ch1 = Chapter;
ch1.Title = 'Preprocessing report';

sec1.Title = 'Bad channel rejection';
append(ch1,sec1);
para = Paragraph('Bad channels were automatically detected with clean_rawdata. By default flat channels, channels with high frequency noise and channels with poor predictability are rejected. In the following plot,  bad channels are marked in blue.');
% Add here figure of bad channels

bc_plot = plot_badchannels_singlestudy(params, bidsID);
bc_f = Figure(bc_plot);
append(ch1,para);
add(ch1, bc_f);


sec2.Title = 'ICA. Independent component classification';
append(ch1,sec2);
[IC_plot, ic_kept, ic_kept_pc] = plot_ICs_singlestudy(params, bidsID);
para = Paragraph('Artefactual independent components were detected automatically with ICLabel. By default, components whose probability of being ''Muscle'' or ''Eye'' are higher than 80% are marked as artifactual and substracted from the data.');
append(ch1,para);
f_ic = Figure(IC_plot);
add(ch1, f_ic);
para = Paragraph([num2str(ic_kept) ' components were kept (' num2str(ic_kept_pc) '%)']);
append(ch1,para);

sec3.Title = 'Bad time segments rejection';
append(ch1,sec3);
[bs_plot, bs_secs, bs_pc] = plot_badtimesegments_singlestudy(params, bidsID);
para = Paragraph(['Bad time segments were detected automatically with ASR. In total ' num2str(bs_pc) ' % of the data (' num2str(bs_secs) ' seconds) were rejected']);

f_bs = Figure(bs_plot);
add(ch1, f_bs);
append(ch1,para);

append(rpt,ch1)
%% EEG features
ch2 = Chapter;
ch2.Title = 'EEG features report';

sec1.Title = 'Power';
append(ch2,sec1);
[power_fig, topoplot_fig] = plot_power(params,bidsID);
f1 = Figure(power_fig);
add(ch2,f1);
f2 = Figure(topoplot_fig);
add(ch2,f2);


sec2.Title = 'Peak frequency';
append(ch2,sec2);
pf_fig = plot_peakfrequency(params,bidsID);
f3 = Figure(pf_fig);
add(ch2,f3);

append(rpt,ch2);
close(rpt);
rptview(rpt);
end