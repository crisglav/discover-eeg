function pdf_report(params,bidsID)

pdf_file = fullfile(params.ReportsPath,[bidsID '.pdf']);
f = figure;
tl = tiledlayout(3,1);
nexttile
bc_plot = plot_badchannels_singlestudy(params, bidsID);


exportgraphics(Figure(bc_plot),pdf_file,"Append",true);

[IC_plot, ic_kept] = plot_ICs_singlestudy(params, bidsID);

add(sec2, Figure(IC_plot));
para = Paragraph([num2str(ic_kept*100, '%.2f') '% IC components were kept.']);
append(sec2,para);

sec3 = Section;
sec3.Title = 'Bad time segments rejection';
append(ch1,sec3);
[bs_plot, bs_length, recording_length] = plot_badtimesegments_singlestudy(params, bidsID);
para = Paragraph(['Bad time segments were detected automatically with ASR. In total ' num2str(bs_length, '%.1f') ' out of ' num2str(recording_length, '%.1f') ' seconds were rejected (' num2str((bs_length/recording_length)*100, '%.2f') ' % of the data).']);
append(sec3,para);
add(sec3, Figure(bs_plot));

append(rpt,ch1)
%% EEG features
ch2 = Chapter;
ch2.Title = 'EEG features report';
sec1 = Section;
sec1.Title = 'Power (sensor space)';
append(ch2,sec1);
para = Paragraph('Power spectrum was computed with dpss multitapers and a frequency smoothing of 1 Hz by default. Average power across all epochs and channels is depicted below. In the legend, the percentage of relative power for each frequency band is included in brackets.');
append(sec1,para);
[power_plot, topo_plot] = plot_power(params,bidsID);
add(sec1,Figure(power_plot));
para = Paragraph('Spatial distribution of power per frequency band. Average power across epochs is depicted for each frequency band.');
append(sec1,para);
add(sec1,Figure(topo_plot));

sec2 = Section;
sec2.Title = 'Peak frequency (sensor space)';
append(ch2,sec2);
para = Paragraph('Peak frequency was computed in the PSD averaged across epochs and channels in the alpha range with two different methods: frequency at the maximum of the peak and the center of gravity of the power spectrum in the alpha range.');
append(sec2,para);
pf_plot = plot_peakfrequency(params,bidsID);
add(sec2,Figure(pf_plot));

sec3 = Section;
sec3.Title = 'Power (source space)';
append(ch2,sec3);
para = Paragraph('Power at the 400 source ROIs was estimated by projecting the frequency-band specific spatial filter to the sensor-space band-pass filtered data. For visualization the 400 ROIs were interpolated to a standard cortical surface.');
append(sec3,para);
power_source_plot = plot_power_source(params,bidsID);
add(sec3,Figure(power_source_plot));

sec4 = Section;
sec4.Title = 'Connectivity (source space)';
append(ch2,sec4);
para = Paragraph('Functional connectivity was estimated on the 400 source ROIs with two measures: debiased weighted PLI (phase-based) and orthogonalised Amplitude Envelope Correlation (amplitude-based).');
append(sec4,para);
para = Paragraph('Connectivity matrices are organized per functional network. The left (up) network submodule corresponds to ROIs of that network located on the left brain hemisphere. The right (down) network submodule corresponds to ROIs of that network located on the right hemisphere.');
append(sec4,para);
dwpli_plot = plot_connectivity(params,bidsID, 'dwpli');
aec_plot = plot_connectivity(params,bidsID, 'aec');
sec41 = Section;
sec41.Title = 'debiased weighted Phase Lag Index';
append(sec4,sec41);
add(sec41,Figure(dwpli_plot));
sec42 = Section;
sec42.Title = 'Amplitude Envelope Correlation';
append(sec4,sec42);
add(sec42,Figure(aec_plot));

sec5 = Section;
sec5.Title = 'Network characterization - graph theory measures based on connectivity';
append(ch2,sec5);
para = Paragraph('For each connectivity measure and frequency band, 2 local graph measures (degree and global clustering coefficient) and 3 global graph measures (global clustering coefficient, global efficiency and smallworldness) were computed.');
append(sec5,para);
para = Paragraph('By default the connectivity matrix is thresholded keeping the 20% of strongest connections and binarized for computing the graph measures.');
append(sec5,para);
sec51 = Section;
sec51.Title = 'dwPLI - Network characterization';
append(sec5,sec51);
[degree_plot, cc_plot, global_plot] = plot_graph_measures(params,bidsID,'dwpli');
add(sec51,Figure(degree_plot));
add(sec51,Figure(cc_plot));
add(sec51,Figure(global_plot));
sec52 = Section;
sec52.Title = 'AEC - Network characterization';
append(sec5,sec52);
[degree_plot, cc_plot, global_plot] = plot_graph_measures(params,bidsID,'aec');
add(sec52,Figure(degree_plot));
add(sec52,Figure(cc_plot));
add(sec52,Figure(global_plot));

append(rpt,ch2);

close(rpt);
% rptview(rpt);
end