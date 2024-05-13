function [data_prop_amp] = calc_wallnormal(x_crit,x_crit_val,refine_max,nsamp,index_harm,data_prop_ampt,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if strcmp(x_crit,'max')

% Max value wall-normal
for k=1:nsamp
    for i=1:numel(index_harm)
        if strcmp (refine_max,'true')
            for j=1:numel(index_span)
                [data_prop_amp_coa(i,j,k),index_xmax(j)]=max(data_prop_ampt(i,:,j,k),[],2);
            
                xmax_fine=linspace(x(index_xmax(j)-3),x(index_xmax(j)+3),numel(x));
                data_prop_amp_fine=spline(x,data_prop_ampt(i,:,j,k),xmax_fine);
                data_prop_amp(i,j,k)=max(data_prop_amp_fine);
            end
        else
            data_prop_amp(i,:,k)=max(data_prop_ampt(i,:,:,k),[],2);
        end
    end
end

elseif strcmp(x_crit,'value')

[~,index_x_crit_val]=min(abs(x_crit_val-x)); 
for k=1:nsamp
    for i=1:numel(index_harm)
        data_prop_amp(i,:,k)=data_prop_ampt(i,index_x_crit_val,:,k);
    end
end

end

end

