function [specy] = load_spectra_new(path_case,var,f,count,nx,index_span,nz,specy)
%LOAD_SPECTRA Summary of this function goes here
%   Detailed explanation goes here

for i=1:numel(index_span)
    fname = sprintf(strcat(path_case,'results/fft/fft_%s_span%i_%07d.bin'), var,index_span(i),f);
    fileID = fopen(fname);
    data2 = fread(fileID,'double');
    fclose(fileID);

    data = reshape(data2,2,nx,nz);

    specy(count,:,index_span(i)+1,:)=complex(data(1,:,:),data(2,:,:));
end

end

