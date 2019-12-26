% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

close all;

%% Parameters: source

%Number of plane waves
Nsrc = 360;

%Plane wave direction (rad)
d_theta_s = 2*pi/Nsrc;
theta_s = (0:Nsrc-1)*d_theta_s;

%Amplitude
amp_src = ones(1,Nsrc);

%% Memory allocation

%SDR
sdrr_reg = zeros(prm.freq_len,1);
sdrr_rnd = zeros(prm.freq_len,1);
sdrr_gso = zeros(prm.freq_len,1);
sdrr_det = zeros(prm.freq_len,1);
sdrr_mi = zeros(prm.freq_len,1);
sdrr_fs = zeros(prm.freq_len,1);
sdrr_eim = zeros(prm.freq_len,1);

%Condition number in dB
cond_db_reg = zeros(prm.freq_len,1);
cond_db_rnd = zeros(prm.freq_len,1);
cond_db_gso = zeros(prm.freq_len,1);
cond_db_det = zeros(prm.freq_len,1);
cond_db_mi = zeros(prm.freq_len,1);
cond_db_fs = zeros(prm.freq_len,1);
cond_db_eim = zeros(prm.freq_len,1);

for kk=1:prm.freq_len
    fprintf('freq: %d\n',prm.freq(kk));
    
    %% Calculate transfer functions

    H_l2m_reg = H_lc2d(m_reg_idx{kk},l_reg_idx{kk},kk);
    H_l2m_rnd = H_lc2d(m_rnd_idx{kk},l_rnd_idx{kk},kk);
    H_l2m_gso = H_lc2d(m_gso_idx{kk},l_gso_idx{kk},kk);
    H_l2m_det = H_lc2d(m_det_idx{kk},l_det_idx{kk},kk);
    H_l2m_mi = H_lc2d(m_mi_idx{kk},l_mi_idx{kk},kk);
    H_l2m_fs = H_lc2d(m_fs_idx{kk},l_fs_idx{kk},kk);
    H_l2m_eim = H_lc2d(m_eim_idx{kk},l_eim_idx{kk},kk);

    %% Desired sound field 

    des = zeros(des_prm.num,Nsrc);
    des_m_reg = zeros(m_reg_num(kk),Nsrc);
    des_m_rnd = zeros(m_rnd_num(kk),Nsrc);
    des_m_gso = zeros(m_gso_num(kk),Nsrc);
    des_m_det = zeros(m_det_num(kk),Nsrc);
    des_m_mi = zeros(m_mi_num(kk),Nsrc);
    des_m_fs = zeros(m_fs_num(kk),Nsrc);
    des_m_eim = zeros(m_eim_num(kk),Nsrc);

    for is=1:Nsrc
        %Target
        H_s2d = amp_src(is).*exp(1i*prm.k(kk)*(des_prm.pos(:,1)*cos(theta_s(is))+des_prm.pos(:,2)*sin(theta_s(is))));
        des(:,is) = des(:,is) + H_s2d;

        %Reg
        H_s2m_reg = amp_src(is).*exp(1i*prm.k(kk)*(m_reg_pos{kk}(:,1)*cos(theta_s(is))+m_reg_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_reg(:,is) = des_m_reg(:,is) + H_s2m_reg;

        %Rand
        H_s2m_rnd = amp_src(is).*exp(1i*prm.k(kk)*(m_rnd_pos{kk}(:,1)*cos(theta_s(is))+m_rnd_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_rnd(:,is) = des_m_rnd(:,is) + H_s2m_rnd;

        %GSO
        H_s2m_gso = amp_src(is).*exp(1i*prm.k(kk)*(m_gso_pos{kk}(:,1)*cos(theta_s(is))+m_gso_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_gso(:,is) = des_m_gso(:,is) + H_s2m_gso;

        %Det
        H_s2m_det = amp_src(is).*exp(1i*prm.k(kk)*(m_det_pos{kk}(:,1)*cos(theta_s(is))+m_det_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_det(:,is) = des_m_det(:,is) + H_s2m_det;

        %MI
        H_s2m_mi = amp_src(is).*exp(1i*prm.k(kk)*(m_mi_pos{kk}(:,1)*cos(theta_s(is))+m_mi_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_mi(:,is) = des_m_mi(:,is) + H_s2m_mi;

        %FS
        H_s2m_fs = amp_src(is).*exp(1i*prm.k(kk)*(m_fs_pos{kk}(:,1)*cos(theta_s(is))+m_fs_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_fs(:,is) = des_m_fs(:,is) + H_s2m_fs;

        %EIM
        H_s2m_eim = amp_src(is).*exp(1i*prm.k(kk)*(m_eim_pos{kk}(:,1)*cos(theta_s(is))+m_eim_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_eim(:,is) = des_m_eim(:,is) + H_s2m_eim;
    end

    %% Pressure matching

    %Reg
    drv_reg = pinv(H_l2m_reg)*des_m_reg;
    cond_db_reg(kk) = 20*log10(cond(H_l2m_reg));

    %Rand
    drv_rnd = pinv(H_l2m_rnd)*des_m_rnd;
    cond_db_rnd(kk) = 20*log10(cond(H_l2m_rnd));

    %GSO
    drv_gso = pinv(H_l2m_gso)*des_m_gso;
    cond_db_gso(kk) = 20*log10(cond(H_l2m_gso));

    %Det
    drv_det = pinv(H_l2m_det)*des_m_det;
    cond_db_det(kk) = 20*log10(cond(H_l2m_det));

    %MI
    drv_mi = pinv(H_l2m_mi)*des_m_mi;
    cond_db_mi(kk) = 20*log10(cond(H_l2m_mi));

    %FS
    drv_fs = pinv(H_l2m_fs)*des_m_fs;
    cond_db_fs(kk) = 20*log10(cond(H_l2m_fs));

    %EIM
    drv_eim = pinv(H_l2m_eim)*des_m_eim;
    cond_db_eim(kk) = 20*log10(cond(H_l2m_eim));

    %% Reproduction

    fprintf('Reproduction...\n');

    amp_src_mat = ones(rep_prm.num,1)*amp_src;
    ideal = amp_src_mat.*exp(1i*prm.k(kk)*(rep_prm.pos(:,1)*cos(theta_s)+rep_prm.pos(:,2)*sin(theta_s)));

    H_l2r_reg = H_lc2r(:,l_reg_idx{kk},kk);
    H_l2r_rnd = H_lc2r(:,l_rnd_idx{kk},kk);
    H_l2r_gso = H_lc2r(:,l_gso_idx{kk},kk);
    H_l2r_det = H_lc2r(:,l_det_idx{kk},kk);
    H_l2r_mi = H_lc2r(:,l_mi_idx{kk},kk);
    H_l2r_fs = H_lc2r(:,l_fs_idx{kk},kk);
    H_l2r_eim = H_lc2r(:,l_eim_idx{kk},kk);

    %Synthesized pressure distribution
    syn_reg = H_l2r_reg*drv_reg;
    syn_rnd = H_l2r_rnd*drv_rnd;
    syn_gso = H_l2r_gso*drv_gso;
    syn_det = H_l2r_det*drv_det;
    syn_mi = H_l2r_mi*drv_mi;
    syn_fs = H_l2r_fs*drv_fs;
    syn_eim = H_l2r_eim*drv_eim;

    %% Evaluation

    %SDR
    sdrr_reg_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_reg(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_rnd_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_rnd(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_gso_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_gso(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_det_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_det(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_mi_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_mi(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_fs_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_fs(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_eim_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_eim(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));

    sdrr_reg(kk) = mean(sdrr_reg_mat);
    sdrr_rnd(kk) = mean(sdrr_rnd_mat);
    sdrr_gso(kk) = mean(sdrr_gso_mat);
    sdrr_det(kk) = mean(sdrr_det_mat);
    sdrr_mi(kk) = mean(sdrr_mi_mat);
    sdrr_fs(kk) = mean(sdrr_fs_mat);
    sdrr_eim(kk) = mean(sdrr_eim_mat);

    fprintf('[SDR] Reg: %f, Rand: %f, GSO: %f, Det: %f, MI: %f, FS: %f, EIM: %f\n',sdrr_reg(kk),sdrr_rnd(kk),sdrr_gso(kk),sdrr_det(kk),sdrr_mi(kk),sdrr_fs(kk),sdrr_eim(kk));
    fprintf('[Cond num] Reg: %f, Rand: %f, GSO: %f, Det: %f, MI: %f, FS: %f, EIM: %f\n',cond_db_reg(kk),cond_db_rnd(kk),cond_db_gso(kk),cond_db_mi(kk),cond_db_det(kk),cond_db_fs(kk),cond_db_eim(kk));

    %% Draw figures

    if ismember(prm.freq(kk), figsave_freq) == 1
        idx_plt = 220; %Index for plot

        ideal_plt = reshape(ideal(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_reg_plt = reshape(syn_reg(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_rnd_plt = reshape(syn_rnd(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_gso_plt = reshape(syn_gso(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_det_plt = reshape(syn_det(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_mi_plt = reshape(syn_mi(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_fs_plt = reshape(syn_fs(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_eim_plt = reshape(syn_eim(:,idx_plt),rep_prm.Nx,rep_prm.Ny);

        err_reg_plt = 10*log10((abs(syn_reg_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_rnd_plt = 10*log10((abs(syn_rnd_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_gso_plt = 10*log10((abs(syn_gso_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_det_plt = 10*log10((abs(syn_det_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_mi_plt = 10*log10((abs(syn_mi_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_fs_plt = 10*log10((abs(syn_fs_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_eim_plt = 10*log10((abs(syn_eim_plt-ideal_plt).^2)./(abs(ideal_plt).^2));

        sdrr_reg_plt = sdrr_reg_mat(idx_plt);
        sdrr_rnd_plt = sdrr_rnd_mat(idx_plt);
        sdrr_gso_plt = sdrr_gso_mat(idx_plt);
        sdrr_det_plt = sdrr_det_mat(idx_plt);
        sdrr_mi_plt = sdrr_mi_mat(idx_plt);
        sdrr_fs_plt = sdrr_fs_mat(idx_plt);
        sdrr_eim_plt = sdrr_eim_mat(idx_plt);

        %Synthesized pressure distribution
        fig_pos = [0, 0, 350, 280];
        z_range = [-1.5, 1.5];
        xx = [min(rep_prm.pos(:,1)),max(rep_prm.pos(:,1))];
        yy = [min(rep_prm.pos(:,2)),max(rep_prm.pos(:,2))];

        %Target
        fig(11)=figure(11);
        set(fig(11),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(ideal_plt).');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %Reg
        fig(12)=figure(12);
        set(fig(12),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_reg_plt).');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %Rand
        fig(13)=figure(13);
        set(fig(13),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_rnd_plt).');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %GSO
        fig(14)=figure(14);
        set(fig(14),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_gso_plt).');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %Det
        fig(15)=figure(15);
        set(fig(15),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_det_plt).');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %MI
        fig(16)=figure(16);
        set(fig(16),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_mi_plt).');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %FS
        fig(17)=figure(17);
        set(fig(17),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_fs_plt).');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %EIM
        fig(18)=figure(18);
        set(fig(18),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_eim_plt).');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %Normalized error distribution
        z_range = [-40, 10];

        %Reg
        fig(21)=figure(21);
        set(fig(21),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_reg_plt.');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %Rand
        fig(22)=figure(22);
        set(fig(22),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_rnd_plt.');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %GSO
        fig(23)=figure(23);
        set(fig(23),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_gso_plt.');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %Det
        fig(24)=figure(24);
        set(fig(24),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_det_plt.');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %MI
        fig(25)=figure(25);
        set(fig(25),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_mi_plt.');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %FS
        fig(26)=figure(26);
        set(fig(26),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_fs_plt.');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %EIM
        fig(27)=figure(27);
        set(fig(27),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_eim_plt.');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %Loudspeaker output power
        fig_pos = [0, 0, 600, 250];

        fig(31)=figure(31);
        set(fig(31),'Position',fig_pos);
        boxplot([10*log10(abs(drv_reg(:)).^2),10*log10(abs(drv_rnd(:)).^2),10*log10(abs(drv_gso(:)).^2),10*log10(abs(drv_det(:)).^2),10*log10(abs(drv_mi(:)).^2),10*log10(abs(drv_fs(:)).^2),10*log10(abs(drv_eim(:)).^2)],'Labels',{'Reg','Rand','GSO','Det','MI','FS','EIM'},'Whisker',1.5);
        ylabel('Output power (dB)');

        drawnow;

        if figsave_flg==1
            fname = cell(15,1);

            fname{1} = sprintf('%s/fig7a',dir_results);
            fname{2} = sprintf('%s/fig7b',dir_results);
            fname{3} = sprintf('%s/fig7c',dir_results);
            fname{4} = sprintf('%s/fig7d',dir_results);
            fname{5} = sprintf('%s/fig7e',dir_results);
            fname{6} = sprintf('%s/fig7f',dir_results);
            fname{7} = sprintf('%s/fig7g',dir_results);
            fname{8} = sprintf('%s/fig8a',dir_results);
            fname{9} = sprintf('%s/fig8b',dir_results);
            fname{10} = sprintf('%s/fig8c',dir_results);
            fname{11} = sprintf('%s/fig8d',dir_results);
            fname{12} = sprintf('%s/fig8e',dir_results);
            fname{13} = sprintf('%s/fig8f',dir_results);
            fname{14} = sprintf('%s/fig8g',dir_results);
            fname{15} = sprintf('%s/fig9',dir_results);
            save_figures([fig(12:18),fig(21:27),fig(31)], fname, fig_type);
        end
    end
end

%% Draw figures
    
fig_pos = [0, 0, 500, 400];

%Number of sources/sensors
fig(1)=figure(1);
set(fig(1),'Position',fig_pos);
plot(prm.freq,m_eim_num,'-ok');
grid on;
xlim([20,1600]);
ylim([0,50]);
xlabel('Frequency (Hz)'); ylabel('Number of sources/sensors');

fig_pos = [0, 0, 800, 350];

%SDR
fig(2)=figure(2);
set(fig(2),'Position',fig_pos);
hold on;
plot(prm.freq,sdrr_reg,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,sdrr_rnd,'--x','Color',red,'MarkerSize',6);
plot(prm.freq,sdrr_gso,'-.*','Color',yellow,'MarkerSize',6);
plot(prm.freq,sdrr_det,'-^','Color',purple,'MarkerSize',6);
hold off;
grid on;
xlim([20,1600]);
ylim([0,60]);
xlabel('Frequency (Hz)'); ylabel('SDR (dB)');
legend('Reg','Rand','GSO','Det','Location','southwest');

fig(3)=figure(3);
set(fig(3),'Position',fig_pos);
hold on;
plot(prm.freq,sdrr_reg,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,sdrr_mi,'--v','Color',green,'MarkerSize',6);
plot(prm.freq,sdrr_fs,'-.s','Color',sky,'MarkerSize',6);
plot(prm.freq,sdrr_eim,'-+','Color',brown,'MarkerSize',6);
hold off;
grid on;
xlim([20,1600]);
ylim([0,60]);
xlabel('Frequency (Hz)'); ylabel('SDR (dB)');
legend('Reg','MI','FS','EIM','Location','southwest');

%Condition number 
fig(4)=figure(4);
set(fig(4),'Position',fig_pos);
hold on;
plot(prm.freq,cond_db_reg,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,cond_db_rnd,'--x','Color',red,'MarkerSize',6);
plot(prm.freq,cond_db_gso,'-.*','Color',yellow,'MarkerSize',6);
plot(prm.freq,cond_db_det,'-^','Color',purple,'MarkerSize',6);
hold off;
grid on;
xlim([20,1600]);
ylim([40,150]);
xlabel('Frequency (Hz)'); ylabel('Condition number (dB)');
legend('Reg','Rand','GSO','Det','Location','northeast');

fig(5)=figure(5);
set(fig(5),'Position',fig_pos);
hold on;
plot(prm.freq,cond_db_reg,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,cond_db_mi,'--v','Color',green,'MarkerSize',6);
plot(prm.freq,cond_db_fs,'-.s','Color',sky,'MarkerSize',6);
plot(prm.freq,cond_db_eim,'-+','Color',brown,'MarkerSize',6);
hold off;
grid on;
xlim([20,1600]);
ylim([40,150]);
xlabel('Frequency (Hz)'); ylabel('Condition number (dB)');
legend('Reg','MI','FS','EIM','Location','northeast');

if figsave_flg==1
    fname = cell(5,1);

    fname{1} = sprintf('%s/fig3',dir_results);
    fname{2} = sprintf('%s/fig4a',dir_results);
    fname{3} = sprintf('%s/fig4b',dir_results);
    fname{4} = sprintf('%s/fig5a',dir_results);
    fname{5} = sprintf('%s/fig5b',dir_results);
    save_figures(fig(1:5), fname, fig_type);
end
