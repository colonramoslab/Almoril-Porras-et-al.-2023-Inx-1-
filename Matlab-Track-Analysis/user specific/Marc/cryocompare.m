memory 
tic
cryo = ExperimentSet.fromMatFiles(fullfile('\\labnas1\Share\David\Extracted\Spatial\N2\18-23GradientC15\OutputFiles\', 'matfiles','cryo_spatial'));
toc
memory
tic
ctemp = ExperimentSet.fromMatFiles(fullfile('\\labnas1\Share\David\Extracted\Temporal\N2G15\Good\all\', 'matfiles','cryo_temporal'));
toc
memory
hcr1 = cryo.makeHistogram('covRatio', 1:0.1:4);
hcr2 = ctemp.makeHistogram('covRatio', 1:0.1:4);

plot (1:0.1:4, hcr1/sum(hcr1), 1:0.1:4, hcr2/sum(hcr2));