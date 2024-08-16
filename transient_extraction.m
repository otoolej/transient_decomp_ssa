%to open and read the pre-processed NIRS; this file might be mydata1.txt or
%mydata2.txt 
data = readtable('mydata.txt', 'FileType', 'text');
data = table2array(data)
%The first column is time in seconds and the second one is preprocessed
%NIRS
sat = data(:,2)
%to fill in the NaN values in sat
sat = fillmissing(sat, 'linear')
%to decompose the transients; developed by O'Toole JM. Dempsey EM, Boylan 
%GB (2018) 'Extracting transients from cerebral oxygenation signals of 
%preterm infants: a new singular-spectrum analysis method' in Int Conf 
%IEEE Eng Med Biol Society (EMBC), IEEE, pp. 5882--5885 DOI:10.1109/EMBC.2018.8513523
params = decomp_PARAMS
fs = 1 / 6;
db_plot = false;
y = shorttime_iter_SSA_decomp(sat, fs, params, db_plot);
transient = y.component.'
nirs = y.nirs.'
without_component = nirs - transient + nanmean(sat)
%save the different components as text files
writematrix(without_component, 'filtered.txt', 'Delimiter', 'tab')
writematrix(transient, 'transient.txt', 'Delimiter', 'tab')
