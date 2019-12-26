% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

close all;

%Regularization parameters
prm.reg_num = 100;
prm.reg_min = -6.0;
prm.reg_max = 1.0;
prm.reg = 10.^(((1:prm.reg_num)-1)/prm.reg_num*(prm.reg_max-prm.reg_min)+prm.reg_min);

%% Parameters: source

%Number of plane waves
Nsrc = 360;

%Plane wave direction (rad)
d_theta_s = 2*pi/Nsrc;
theta_s = (0:Nsrc-1)*d_theta_s;

%Amplitude
amp_src = ones(1,Nsrc);

%% Position error 

%Standard deviation
l_pos_err_std = 0.01;
m_pos_err_std = 0.01;

%Number of trials
prm.num_tr = 10;

if rev_prm.flg==0
    %Freefield
    l_pos_err = cell(prm.num_tr,1);
    m_pos_err = cell(prm.num_tr,1);
    
    pos_rd_e = cell(prm.num_tr,1);
    for ii=1:prm.num_tr
        l_pos_err{ii} = l_pos_err_std*randn(lc_prm.num,2);
        m_pos_err{ii} = m_pos_err_std*randn(des_prm.num,2); 
    end
    
    for ii=1:prm.num_tr
        lc_prm.pos_e_cell{ii} = lc_prm.pos + l_pos_err{ii};
        des_prm.pos_e_cell{ii} = des_prm.pos + m_pos_err{ii};

        pos_rd_e{ii} = [rep_prm.pos; des_prm.pos_e_cell{ii}];
    end
    
    pos_rd_e_num = size(pos_rd_e{1},1);
    
    H_lc2rd_e_cell = cell(prm.num_tr,1);
    for ii=1:prm.num_tr
        H_lc2rd_e_cell{ii} = zeros(pos_rd_e_num, lc_prm.num, prm.freq_len);
        for kk=1:prm.freq_len
            H_lc2rd_e_cell{ii}(:,:,kk) = green2d(repmat(pos_rd_e{ii}(:,1),1,lc_prm.num),repmat(pos_rd_e{ii}(:,2),1,lc_prm.num),repmat(lc_prm.pos_e_cell{ii}(:,1).',pos_rd_e_num,1),repmat(lc_prm.pos_e_cell{ii}(:,2).',pos_rd_e_num,1),prm.k(kk));
        end
    end
else
    %Reverberant (load files)
    H_lc2rd_e_cell = cell(prm.num_tr,1);
    pos_rd_e_num = size([rep_prm.pos; des_prm.pos],1);
    for ii=1:prm.num_tr
        H_lc2rd_e_cell{ii} = zeros(pos_rd_e_num, lc_prm.num, prm.freq_len);
        for kk=1:prm.freq_len
            fname_imp_lc2rd_e_f = sprintf('../data/imp_lc2rd_e_abs%.2f_l%d_%dx%d_std%.2f_r%d_%dx%d_d%d_%dx%d_std%.2f/imp_lc2rd_e%02d_f%d.mat',rev_prm.abs_ratio,lc_prm.num,lc_prm.len_x*100,lc_prm.len_y*100,l_pos_err_std,rep_prm.num,rep_prm.len_x*100,rep_prm.len_y*100,des_prm.num,des_prm.len_x*100,des_prm.len_y*100,m_pos_err_std,ii,prm.freq(kk));
            load(fname_imp_lc2rd_e_f,'lc_prm','des_prm','H_lc2rd_e_f');
            H_lc2rd_e_cell{ii}(:,:,kk) = H_lc2rd_e_f;
        end
    end
end

%% Memory allocation

%Regularization parameter
regprm_reg = zeros(prm.freq_len,prm.num_tr);
regprm_rnd = zeros(prm.freq_len,prm.num_tr);
regprm_gso = zeros(prm.freq_len,prm.num_tr);
regprm_det = zeros(prm.freq_len,prm.num_tr);
regprm_mi = zeros(prm.freq_len,prm.num_tr);
regprm_fs = zeros(prm.freq_len,prm.num_tr);
regprm_eim = zeros(prm.freq_len,prm.num_tr);

regprm_reg_idx = zeros(prm.freq_len,prm.num_tr);
regprm_rnd_idx = zeros(prm.freq_len,prm.num_tr);
regprm_gso_idx = zeros(prm.freq_len,prm.num_tr);
regprm_det_idx = zeros(prm.freq_len,prm.num_tr);
regprm_mi_idx = zeros(prm.freq_len,prm.num_tr);
regprm_fs_idx = zeros(prm.freq_len,prm.num_tr);
regprm_eim_idx = zeros(prm.freq_len,prm.num_tr);

%SDR
sdrr_reg_opt = zeros(prm.freq_len,prm.num_tr);
sdrr_rnd_opt = zeros(prm.freq_len,prm.num_tr);
sdrr_gso_opt = zeros(prm.freq_len,prm.num_tr);
sdrr_det_opt = zeros(prm.freq_len,prm.num_tr);
sdrr_mi_opt = zeros(prm.freq_len,prm.num_tr);
sdrr_fs_opt = zeros(prm.freq_len,prm.num_tr);
sdrr_eim_opt = zeros(prm.freq_len,prm.num_tr);

%Condition number in dB
cond_db_reg_opt = zeros(prm.freq_len,prm.num_tr);
cond_db_rnd_opt = zeros(prm.freq_len,prm.num_tr);
cond_db_gso_opt = zeros(prm.freq_len,prm.num_tr);
cond_db_det_opt = zeros(prm.freq_len,prm.num_tr);
cond_db_mi_opt = zeros(prm.freq_len,prm.num_tr);
cond_db_fs_opt = zeros(prm.freq_len,prm.num_tr);
cond_db_eim_opt = zeros(prm.freq_len,prm.num_tr);

for ii=1:prm.num_tr
    fprintf('Trial: %d\n',ii);
    
    H_lc2r_e = H_lc2rd_e_cell{ii}(1:rep_prm.num,:,:);
    H_lc2d_e = H_lc2rd_e_cell{ii}(rep_prm.num+1:rep_prm.num+des_prm.num,:,:);

    cond_db_reg = zeros(prm.reg_num,prm.freq_len);
    cond_db_rnd = zeros(prm.reg_num,prm.freq_len);
    cond_db_gso = zeros(prm.reg_num,prm.freq_len);
    cond_db_det = zeros(prm.reg_num,prm.freq_len);
    cond_db_mi = zeros(prm.reg_num,prm.freq_len);
    cond_db_fs = zeros(prm.reg_num,prm.freq_len);
    cond_db_eim = zeros(prm.reg_num,prm.freq_len);

    for kk=1:prm.freq_len
        fprintf('freq: %d\n',prm.freq(kk));

        %% Calculate transfer functions

        H_l2m_reg = H_lc2d_e(m_reg_idx,l_reg_idx,kk);
        H_l2m_rnd = H_lc2d_e(m_rnd_idx,l_rnd_idx,kk);
        H_l2m_gso = H_lc2d_e(m_gso_idx,l_gso_idx,kk);
        H_l2m_det = H_lc2d_e(m_det_idx,l_det_idx,kk);
        H_l2m_mi = H_lc2d_e(m_mi_idx,l_mi_idx,kk);
        H_l2m_fs = H_lc2d_e(m_fs_idx,l_fs_idx,kk);
        H_l2m_eim = H_lc2d_e(m_eim_idx,l_eim_idx,kk);

        %% Desired sound field 

        des = zeros(des_prm.num,Nsrc);
        des_m_reg = zeros(m_reg_num,Nsrc);
        des_m_rnd = zeros(m_rnd_num,Nsrc);
        des_m_gso = zeros(m_gso_num,Nsrc);
        des_m_det = zeros(m_det_num,Nsrc);
        des_m_mi = zeros(m_mi_num,Nsrc);
        des_m_fs = zeros(m_fs_num,Nsrc);
        des_m_eim = zeros(m_eim_num,Nsrc);

        for is=1:Nsrc
            %Target
            H_s2d = amp_src(is).*exp(1i*prm.k(kk)*(des_prm.pos(:,1)*cos(theta_s(is))+des_prm.pos(:,2)*sin(theta_s(is))));
            des(:,is) = des(:,is) + H_s2d;

            %Reg+EIMi
            H_s2m_reg = amp_src(is).*exp(1i*prm.k(kk)*(m_reg_pos(:,1)*cos(theta_s(is))+m_reg_pos(:,2)*sin(theta_s(is))));
            des_m_reg(:,is) = des_m_reg(:,is) + H_s2m_reg;

            %Rand
            H_s2m_rnd = amp_src(is).*exp(1i*prm.k(kk)*(m_rnd_pos(:,1)*cos(theta_s(is))+m_rnd_pos(:,2)*sin(theta_s(is))));
            des_m_rnd(:,is) = des_m_rnd(:,is) + H_s2m_rnd;

            %GSO
            H_s2m_gso = amp_src(is).*exp(1i*prm.k(kk)*(m_gso_pos(:,1)*cos(theta_s(is))+m_gso_pos(:,2)*sin(theta_s(is))));
            des_m_gso(:,is) = des_m_gso(:,is) + H_s2m_gso;

            %Det
            H_s2m_det = amp_src(is).*exp(1i*prm.k(kk)*(m_det_pos(:,1)*cos(theta_s(is))+m_det_pos(:,2)*sin(theta_s(is))));
            des_m_det(:,is) = des_m_det(:,is) + H_s2m_det;

            %MI
            H_s2m_mi = amp_src(is).*exp(1i*prm.k(kk)*(m_mi_pos(:,1)*cos(theta_s(is))+m_mi_pos(:,2)*sin(theta_s(is))));
            des_m_mi(:,is) = des_m_mi(:,is) + H_s2m_mi;

            %FS
            H_s2m_fs = amp_src(is).*exp(1i*prm.k(kk)*(m_fs_pos(:,1)*cos(theta_s(is))+m_fs_pos(:,2)*sin(theta_s(is))));
            des_m_fs(:,is) = des_m_fs(:,is) + H_s2m_fs;

            %EIM
            H_s2m_eim = amp_src(is).*exp(1i*prm.k(kk)*(m_eim_pos(:,1)*cos(theta_s(is))+m_eim_pos(:,2)*sin(theta_s(is))));
            des_m_eim(:,is) = des_m_eim(:,is) + H_s2m_eim;
        end

        %% Regularization parameter search

        %Loudspeaker driving signals
        drv_reg = zeros(l_reg_num,Nsrc,prm.reg_num);
        drv_rnd = zeros(l_rnd_num,Nsrc,prm.reg_num);
        drv_gso = zeros(l_gso_num,Nsrc,prm.reg_num);
        drv_det = zeros(l_det_num,Nsrc,prm.reg_num);
        drv_mi = zeros(l_mi_num,Nsrc,prm.reg_num);
        drv_fs = zeros(l_fs_num,Nsrc,prm.reg_num);
        drv_eim = zeros(l_eim_num,Nsrc,prm.reg_num);

        %Synthsized pressure distribution
        syn_reg = zeros(rep_prm.num,Nsrc,prm.reg_num);
        syn_rnd = zeros(rep_prm.num,Nsrc,prm.reg_num);
        syn_gso = zeros(rep_prm.num,Nsrc,prm.reg_num);
        syn_det = zeros(rep_prm.num,Nsrc,prm.reg_num);
        syn_mi = zeros(rep_prm.num,Nsrc,prm.reg_num);
        syn_fs = zeros(rep_prm.num,Nsrc,prm.reg_num);
        syn_eim = zeros(rep_prm.num,Nsrc,prm.reg_num);

        fprintf('Reproduction...\n');

        %SDR
        sdrr_reg_mat = zeros(prm.reg_num,Nsrc);
        sdrr_rnd_mat = zeros(prm.reg_num,Nsrc);
        sdrr_gso_mat = zeros(prm.reg_num,Nsrc);
        sdrr_det_mat = zeros(prm.reg_num,Nsrc);
        sdrr_mi_mat = zeros(prm.reg_num,Nsrc);
        sdrr_fs_mat = zeros(prm.reg_num,Nsrc);
        sdrr_eim_mat = zeros(prm.reg_num,Nsrc);

        for rr=1:prm.reg_num
            %% Pressure matching

            lambda = prm.reg(rr);

            %Reg+EIMi
            drv_reg(:,:,rr) = pinv(H_l2m_reg'*H_l2m_reg+lambda*eye(l_reg_num))*H_l2m_reg'*des_m_reg;
            cond_db_reg(rr,kk) = 10*log10(cond(H_l2m_reg'*H_l2m_reg+lambda*eye(l_reg_num)));

            %Rand
            drv_rnd(:,:,rr) = pinv(H_l2m_rnd'*H_l2m_rnd+lambda*eye(l_rnd_num))*H_l2m_rnd'*des_m_rnd;
            cond_db_rnd(rr,kk) = 10*log10(cond(H_l2m_rnd'*H_l2m_rnd+lambda*eye(l_rnd_num)));

            %GSO
            drv_gso(:,:,rr) = pinv(H_l2m_gso'*H_l2m_gso+lambda*eye(l_gso_num))*H_l2m_gso'*des_m_gso;
            cond_db_gso(rr,kk) = 10*log10(cond(H_l2m_gso'*H_l2m_gso+lambda*eye(l_gso_num)));

            %Det
            drv_det(:,:,rr) = pinv(H_l2m_det'*H_l2m_det+lambda*eye(l_det_num))*H_l2m_det'*des_m_det;
            cond_db_det(rr,kk) = 10*log10(cond(H_l2m_det'*H_l2m_det+lambda*eye(l_det_num)));

            %MI
            drv_mi(:,:,rr) = pinv(H_l2m_mi'*H_l2m_mi+lambda*eye(l_mi_num))*H_l2m_mi'*des_m_mi;
            cond_db_mi(rr,kk) = 10*log10(cond(H_l2m_mi*H_l2m_mi'+lambda*eye(l_mi_num)));

            %FS
            drv_fs(:,:,rr) = pinv(H_l2m_fs'*H_l2m_fs+lambda*eye(l_fs_num))*H_l2m_fs'*des_m_fs;
            cond_db_fs(rr,kk) = 10*log10(cond(H_l2m_fs'*H_l2m_fs+lambda*eye(l_fs_num)));

            %EIM
            drv_eim(:,:,rr) = pinv(H_l2m_eim'*H_l2m_eim+lambda*eye(l_eim_num))*H_l2m_eim'*des_m_eim;
            cond_db_eim(rr,kk) = 10*log10(cond(H_l2m_eim'*H_l2m_eim+lambda*eye(l_eim_num)));

            %% Reproduction

            amp_src_mat = ones(rep_prm.num,1)*amp_src;
            ideal = amp_src_mat.*exp(1i*prm.k(kk)*(rep_prm.pos(:,1)*cos(theta_s)+rep_prm.pos(:,2)*sin(theta_s)));

            H_l2r_reg = H_lc2r_e(:,l_reg_idx,kk);
            H_l2r_rnd = H_lc2r_e(:,l_rnd_idx,kk);
            H_l2r_gso = H_lc2r_e(:,l_gso_idx,kk);
            H_l2r_det = H_lc2r_e(:,l_det_idx,kk);
            H_l2r_mi = H_lc2r_e(:,l_mi_idx,kk);
            H_l2r_fs = H_lc2r_e(:,l_fs_idx,kk);
            H_l2r_eim = H_lc2r_e(:,l_eim_idx,kk);

            %Synthesized pressure distribution
            syn_reg(:,:,rr) = H_l2r_reg*drv_reg(:,:,rr);
            syn_rnd(:,:,rr) = H_l2r_rnd*drv_rnd(:,:,rr);
            syn_gso(:,:,rr) = H_l2r_gso*drv_gso(:,:,rr);
            syn_det(:,:,rr) = H_l2r_det*drv_det(:,:,rr);
            syn_mi(:,:,rr) = H_l2r_mi*drv_mi(:,:,rr);
            syn_fs(:,:,rr) = H_l2r_fs*drv_fs(:,:,rr);
            syn_eim(:,:,rr) = H_l2r_eim*drv_eim(:,:,rr);

            %% Evaluation

            sdrr_reg_mat(rr,:) = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_reg(rep_prm.idx_des,:,rr)-ideal(rep_prm.idx_des,:)).^2,1));
            sdrr_rnd_mat(rr,:) = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_rnd(rep_prm.idx_des,:,rr)-ideal(rep_prm.idx_des,:)).^2,1));
            sdrr_gso_mat(rr,:) = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_gso(rep_prm.idx_des,:,rr)-ideal(rep_prm.idx_des,:)).^2,1));
            sdrr_det_mat(rr,:) = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_det(rep_prm.idx_des,:,rr)-ideal(rep_prm.idx_des,:)).^2,1));
            sdrr_mi_mat(rr,:) = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_mi(rep_prm.idx_des,:,rr)-ideal(rep_prm.idx_des,:)).^2,1));
            sdrr_fs_mat(rr,:) = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_fs(rep_prm.idx_des,:,rr)-ideal(rep_prm.idx_des,:)).^2,1));
            sdrr_eim_mat(rr,:) = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_eim(rep_prm.idx_des,:,rr)-ideal(rep_prm.idx_des,:)).^2,1));

        end

        %% Evaluation w/ optimal regularization parameter

        %SDR w/ optimal regularization parameter
        [sdrr_reg_opt(kk,ii),regprm_reg_idx(kk,ii)] = max(mean(sdrr_reg_mat,2));
        [sdrr_rnd_opt(kk,ii),regprm_rnd_idx(kk,ii)] = max(mean(sdrr_rnd_mat,2));
        [sdrr_gso_opt(kk,ii),regprm_gso_idx(kk,ii)] = max(mean(sdrr_gso_mat,2));
        [sdrr_det_opt(kk,ii),regprm_det_idx(kk,ii)] = max(mean(sdrr_det_mat,2));
        [sdrr_mi_opt(kk,ii),regprm_mi_idx(kk,ii)] = max(mean(sdrr_mi_mat,2));
        [sdrr_fs_opt(kk,ii),regprm_fs_idx(kk,ii)] = max(mean(sdrr_fs_mat,2));
        [sdrr_eim_opt(kk,ii),regprm_eim_idx(kk,ii)] = max(mean(sdrr_eim_mat,2));

        %Optimal regularization parameter
        regprm_reg(kk,ii) = prm.reg(regprm_reg_idx(kk));
        regprm_rnd(kk,ii) = prm.reg(regprm_rnd_idx(kk));
        regprm_gso(kk,ii) = prm.reg(regprm_gso_idx(kk));
        regprm_det(kk,ii) = prm.reg(regprm_det_idx(kk));
        regprm_mi(kk,ii) = prm.reg(regprm_mi_idx(kk));
        regprm_fs(kk,ii) = prm.reg(regprm_fs_idx(kk));
        regprm_eim(kk,ii) = prm.reg(regprm_eim_idx(kk));

        %Condition number w/ optimal regularization parameter
        cond_db_reg_opt(kk,ii) = cond_db_reg(regprm_reg_idx(kk),kk);
        cond_db_rnd_opt(kk,ii) = cond_db_rnd(regprm_rnd_idx(kk),kk);
        cond_db_gso_opt(kk,ii) = cond_db_gso(regprm_gso_idx(kk),kk);
        cond_db_det_opt(kk,ii) = cond_db_det(regprm_det_idx(kk),kk);
        cond_db_mi_opt(kk,ii) = cond_db_mi(regprm_mi_idx(kk),kk);
        cond_db_fs_opt(kk,ii) = cond_db_fs(regprm_fs_idx(kk),kk);
        cond_db_eim_opt(kk,ii) = cond_db_eim(regprm_eim_idx(kk),kk);

        fprintf('[Reg param] Reg+EIMi: %f, Rand: %f, GSO: %f, Det: %f, MI: %f, FS: %f, EIM: %f\n',regprm_reg(kk,ii),regprm_rnd(kk,ii),regprm_gso(kk,ii),regprm_det(kk,ii),regprm_mi(kk,ii),regprm_fs(kk,ii),regprm_eim(kk,ii))
        fprintf('[SDR] Reg+EIMi: %f, Rand: %f, GSO: %f, Det: %f, MI: %f, FS: %f, EIM: %f\n',sdrr_reg_opt(kk,ii),sdrr_rnd_opt(kk,ii),sdrr_gso_opt(kk,ii),sdrr_det_opt(kk,ii),sdrr_mi_opt(kk,ii),sdrr_fs_opt(kk,ii),sdrr_eim_opt(kk,ii));
        fprintf('[Cond num] Reg+EIMi: %f, Rand: %f, GSO: %f, Det: %f, MI: %f, FS: %f, EIM: %f\n',cond_db_reg_opt(kk,ii),cond_db_rnd_opt(kk,ii),cond_db_gso_opt(kk,ii),cond_db_det_opt(kk,ii),cond_db_mi_opt(kk,ii),cond_db_fs_opt(kk,ii),cond_db_eim_opt(kk,ii));

    end
end
    
%% Draw figures

%Average SDR
sdrr_reg_opt_av = mean(sdrr_reg_opt,2);
sdrr_rnd_opt_av = mean(sdrr_rnd_opt,2);
sdrr_gso_opt_av = mean(sdrr_gso_opt,2);
sdrr_det_opt_av = mean(sdrr_det_opt,2);
sdrr_mi_opt_av = mean(sdrr_mi_opt,2);
sdrr_fs_opt_av = mean(sdrr_fs_opt,2);
sdrr_eim_opt_av = mean(sdrr_eim_opt,2);
    
fig_pos = [0, 0, 800, 350];

%SDR
fig(2)=figure(2);
set(fig(2),'Position',fig_pos);
hold on;
plot(prm.freq,sdrr_reg_opt_av,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,sdrr_rnd_opt_av,'--x','Color',red,'MarkerSize',6);
plot(prm.freq,sdrr_gso_opt_av,'-.*','Color',yellow,'MarkerSize',6);
plot(prm.freq,sdrr_det_opt_av,'-^','Color',purple,'MarkerSize',6);
hold off;
grid on;
xlim([20,1600]);
ylim([0,50]); 
xlabel('Frequency (Hz)'); ylabel('SDR (dB)');
legend('Reg+EIMi','Rand','GSO','Det','Location','northeast');

fig(3)=figure(3);
set(fig(3),'Position',fig_pos);
hold on;
plot(prm.freq,sdrr_reg_opt_av,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,sdrr_mi_opt_av,'--v','Color',green,'MarkerSize',6);
plot(prm.freq,sdrr_fs_opt_av,'-.s','Color',sky,'MarkerSize',6);
plot(prm.freq,sdrr_eim_opt_av,'-+','Color',brown,'MarkerSize',6);
hold off;
grid on;
xlim([20,1600]);
ylim([0,50]); 
xlabel('Frequency (Hz)'); ylabel('SDR (dB)');
legend('Reg+EIMi','MI','FS','EIM','Location','northeast');

%SDR enlarged
fig(4)=figure(4);
set(fig(4),'Position',fig_pos);
hold on;
plot(prm.freq,sdrr_reg_opt_av,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,sdrr_mi_opt_av,'--v','Color',green,'MarkerSize',6);
plot(prm.freq,sdrr_fs_opt_av,'-.s','Color',sky,'MarkerSize',6);
plot(prm.freq,sdrr_eim_opt_av,'-+','Color',brown,'MarkerSize',6);
hold off;
grid on;
xlim([500,1000]);
xlabel('Frequency (Hz)'); ylabel('SDR (dB)');
legend('Reg+EIMi','MI','FS','EIM','Location','northeast');

if figsave_flg==1
    fname = cell(3,1);

    fname{1} = sprintf('%s/fig22a',dir_results);
    fname{2} = sprintf('%s/fig22b',dir_results);
    fname{3} = sprintf('%s/fig22c',dir_results);
    save_figures(fig(2:4), fname, fig_type);
end
