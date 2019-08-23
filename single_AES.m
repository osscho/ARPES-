close all;
clear all;

        file = sprintf('em3_Fe after PES.gph');
        
       fid = fopen(file,'r');                  % reads the NR of lines ...
                n=0; tline = fgetl(fid);
                while ischar(tline)
                    tline = fgetl(fid);
                    n = n+1;
                end; fclose(fid);
                
        data = dlmread(file,' ',[141 1 n-2 7]);     %exclude header plus last line !
        l=1;
        for i=2:2:size(data); 
            energy(l,1) = data(i-1,2);
            AES_intensity(l,1) = data(i,2);
            beam_current(l,1) = data(i,7);
            AES_normalized_intensity(l,1) = AES_intensity(l,1) / beam_current(l,1);
            l=l+1;
        end;
        
        figure();
        plot(energy,AES_normalized_intensity);
        title('')
        
        % ratio Bi/Sb ratio=5.21*(exp(-(x)/0.79))./(1-exp(-(x)/1.09))