function [EXT_DATA] = TDMS_2_struct( DATA )
% Extract relevant Data from DATA.Data.MeasuredData(k)
% each entry k describes a column of the table, total 12 entries
% Col 1: root (not relavant)    Col 2: some info (not relevant)   Col 3: time info
% Col 4: temperature info       Col 5: voltage V_xy                    Col 6: current
% Col 7: resistance             Col 8: magnetic field             Col 9:
% resistance (unknown, NAN)
% Col 10: coil curr (not relevant)       Col 11: 2 voltage  (unknown,NAN)  Col
% 12: current 2 (unknown,NAN)

% Version: 1.0 C. Schmitz 11.11.2016

EXT_DATA.entries=DATA.Data.MeasuredData(3).Total_Samples; %number of entries
EXT_DATA.time=DATA.Data.MeasuredData(3).Data;
EXT_DATA.temp=DATA.Data.MeasuredData(4).Data;
EXT_DATA.V_xy=DATA.Data.MeasuredData(5).Data;
EXT_DATA.I=DATA.Data.MeasuredData(6).Data;
EXT_DATA.R=DATA.Data.MeasuredData(7).Data;
EXT_DATA.B=DATA.Data.MeasuredData(8).Data;
end

