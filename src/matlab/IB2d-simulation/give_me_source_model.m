function[fs]=give_me_source_model(CL,Nc)

% Model for sniffing
load('peg_indices.mat')
t_step = t_step + 1;

if rem(t_step, save_interval) == 0
    CL_save=CL(peg_ind_start:peg_ind_end,:)';
    dlmwrite('CL_data.csv', CL_save, '-append')
end

fs = zeros(size(CL));
for mer=1:Nc
    for pip=1:length(CL)
        if pip >= peg_ind_start && pip <= peg_ind_end
            fs(pip,mer)=-2*CL(pip,mer);
        else
            fs(pip,mer)=-1*CL(pip,mer);
        end
    end
end
save('peg_indices.mat','t_step','-append')