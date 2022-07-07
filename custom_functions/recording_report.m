function recording_report(params,bidsID)

import mlreportgen.report.*
import mlreportgen.dom.*

f = strcat(bidsID, '_report');
rpt = Report(fullfile(params.reports_folder, f), 'pdf');

% Title page
tp = TitlePage;
tp.Title = ['Report of ' bidsID];
append(rpt,tp);

%% Preprocessing report
ch1 = Chapter;
ch1.Title = 'Preprocessing report';

sec1 = Section;
sec1.Title = 'Bad channel rejection';
append(ch1,sec1);
para = Paragraph('Bad channels were automatically detected with clean_rawdata. By default flat channels, channels with high frequency noise and channels with poor predictability are rejected. In the following plot,  bad channels are marked in blue.');
bc_plot = plot_badchannels_singlestudy(params, bidsID);
append(ch1,para);
add(ch1, Figure(bc_plot));

sec2 = Section;
sec2.Title = 'ICA. Independent component classification';
append(ch1,sec2);
para = Paragraph('Artefactual independent components were detected automatically with ICLabel. By default, components whose probability of being ''Muscle'' or ''Eye'' are higher than 80% are marked as artifactual and substracted from the data.');
append(ch1,para);
para = Paragraph(['In the current run, the threshold used was ' num2str(params.IClabel(2,1)*100, '%.2f') '% for Muscle and ' num2str(params.IClabel(3,1)*100, '%.2f') '% for Eye components.']);
append(ch1,para);
[IC_plot, ic_kept] = plot_ICs_singlestudy(params, bidsID);
add(ch1, Figure(IC_plot));
para = Paragraph([num2str(ic_kept*100, '%.2f') '% IC components were kept.']);
append(ch1,para);

sec3 = Section;
sec3.Title = 'Bad time segments rejection';
append(ch1,sec3);
[bs_plot, bs_length, recording_length] = plot_badtimesegments_singlestudy(params, bidsID);
add(ch1, Figure(bs_plot));
para = Paragraph(['Bad time segments were detected automatically with ASR. In total ' num2str(bs_length, '%.2f') ' out of ' num2str(recording_length, '%.2f') ' seconds were rejected (' num2str((bs_length/recording_length)*100, '%.2f') ' % of the data).']);
append(ch1,para);

append(rpt,ch1)
%% EEG features
ch2 = Chapter;
ch2.Title = 'EEG features report';

sec1 = Section;
sec1.Title = 'Power (sensor space)';
append(ch2,sec1);
[power_plot, topo_plot] = plot_power(params,bidsID);
add(ch2,Figure(power_plot));
add(ch2,Figure(topo_plot));

sec2 = Section;
sec2.Title = 'Peak frequency (sensor space)';
append(ch2,sec2);
pf_plot = plot_peakfrequency(params,bidsID);
add(ch2,Figure(pf_plot));

sec3 = Section;
sec3.Title = 'Power (source space)';
append(ch2,sec3);
para = Paragraph('Power estimates at the 400 ROIs were estimated by projecting the spatial filter to the sensor-space band-pass filtered data. For visualization the 400 ROIs were interpotated to a standard cortical surface.');
append(ch2,para);
power_source_plot = plot_power_source(params,bidsID);
add(ch2,Figure(power_source_plot));

sec4 = Section;
sec4.Title = 'Connectivity (source space)';
append(ch2,sec4);
para = Paragraph('Functional connectivity was estimated on the 400 ROIs at the source level with two measures: debiased weighted PLI (phase-based) and orthogonalised Amplitude Envelope Correlation (amplitude-based).');
append(ch2,para);
para = Paragraph('Connectivity matrices are organized per functional network. The left(upper) submodule corresponds to the left ROIs of that network. The right(down) submodule corresponds to the right ROIs of that network.');
append(ch2,para);
dwpli_plot = plot_connectivity(params,bidsID, 'dwpli');
aec_plot = plot_connectivity(params,bidsID, 'aec');
add(ch2,Figure(dwpli_plot));
add(ch2,Figure(aec_plot));

sec5 = Section;
sec5.Title = 'Network characterization - graph theory measures based on connectivity';
append(ch2,sec5);
para = Paragraph('For each connectivity measure and frequency band, 2 local graph measures (degree and global clustering coefficient) and 3 global graph measures (global clustering coefficient, global efficiency and smallworldness) are computed.');
append(ch2,para);
para = Paragraph('dwPLI - Network characterization');
append(ch2,para);
[degree_plot, cc_plot, global_plot] = plot_graph_measures(params,bidsID,'dwpli');
add(ch2,Figure(degree_plot));
add(ch2,Figure(cc_plot));
add(ch2,Figure(global_plot));
para = Paragraph('AEC - Network characterization');
append(ch2,para);
[degree_plot, cc_plot, global_plot] = plot_graph_measures(params,bidsID,'aec');
add(ch2,Figure(degree_plot));
add(ch2,Figure(cc_plot));
add(ch2,Figure(global_plot));

append(rpt,ch2);

close(rpt);
rptview(rpt);
end