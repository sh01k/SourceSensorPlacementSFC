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
sdrr_regang = zeros(prm.freq_len,1);
sdrr_det = zeros(prm.freq_len,1);
sdrr_mi = zeros(prm.freq_len,1);
sdrr_regangdet = zeros(prm.freq_len,1);
sdrr_regangmi = zeros(prm.freq_len,1);
sdrr_eim = zeros(prm.freq_len,1);

%Condition number in dB
cond_db_reg = zeros(prm.freq_len,1);
cond_db_regang = zeros(prm.freq_len,1);
cond_db_det = zeros(prm.freq_len,1);
cond_db_mi = zeros(prm.freq_len,1);
cond_db_regangdet = zeros(prm.freq_len,1);
cond_db_regangmi = zeros(prm.freq_len,1);
cond_db_eim = zeros(prm.freq_len,1);

for kk=1:prm.freq_len
    fprintf('freq: %d\n',prm.freq(kk));
    
    %% Calculate transfer functions

    H_l2m_reg = H_lc2d(m_reg_idx{kk},l_reg_idx{kk},kk);
    H_l2m_regang = H_lc2d(m_reg_idx{kk},l_regang_idx{kk},kk);
    H_l2m_det = H_lc2d(m_det_idx{kk},l_reg_idx{kk},kk);
    H_l2m_mi = H_lc2d(m_mi_idx{kk},l_reg_idx{kk},kk);
    H_l2m_regangdet = H_lc2d(m_regangdet_idx{kk},l_regang_idx{kk},kk);
    H_l2m_regangmi = H_lc2d(m_regangmi_idx{kk},l_regang_idx{kk},kk);
    H_l2m_eim = H_lc2d(m_eim_idx{kk},l_eim_idx{kk},kk);

    %% Desired sound field 

    des = zeros(des_prm.num,Nsrc);
    des_m_reg = zeros(m_reg_num(kk),Nsrc);
    des_m_regang = zeros(m_reg_num(kk),Nsrc);
    des_m_det = zeros(m_det_num(kk),Nsrc);
    des_m_mi = zeros(m_mi_num(kk),Nsrc);
    des_m_regangdet = zeros(m_regangdet_num(kk),Nsrc);
    des_m_regangmi = zeros(m_regangmi_num(kk),Nsrc);
    des_m_eim = zeros(m_eim_num(kk),Nsrc);

    for is=1:Nsrc
        %Target
        H_s2d = amp_src(is).*exp(1i*prm.k(kk)*(des_prm.pos(:,1)*cos(theta_s(is))+des_prm.pos(:,2)*sin(theta_s(is))));
        des(:,is) = des(:,is) + H_s2d;

        %Reg
        H_s2m_reg = amp_src(is).*exp(1i*prm.k(kk)*(m_reg_pos{kk}(:,1)*cos(theta_s(is))+m_reg_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_reg(:,is) = des_m_reg(:,is) + H_s2m_reg;

        %Reg+EqAng
        H_s2m_regang = amp_src(is).*exp(1i*prm.k(kk)*(m_reg_pos{kk}(:,1)*cos(theta_s(is))+m_reg_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_regang(:,is) = des_m_regang(:,is) + H_s2m_regang;

        %Det
        H_s2m_det = amp_src(is).*exp(1i*prm.k(kk)*(m_det_pos{kk}(:,1)*cos(theta_s(is))+m_det_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_det(:,is) = des_m_det(:,is) + H_s2m_det;

        %MI
        H_s2m_mi = amp_src(is).*exp(1i*prm.k(kk)*(m_mi_pos{kk}(:,1)*cos(theta_s(is))+m_mi_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_mi(:,is) = des_m_mi(:,is) + H_s2m_mi;

        %Det+EqAng
        H_s2m_regangdet = amp_src(is).*exp(1i*prm.k(kk)*(m_regangdet_pos{kk}(:,1)*cos(theta_s(is))+m_regangdet_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_regangdet(:,is) = des_m_regangdet(:,is) + H_s2m_regangdet;

        %MI+EqAng
        H_s2m_regangmi = amp_src(is).*exp(1i*prm.k(kk)*(m_regangmi_pos{kk}(:,1)*cos(theta_s(is))+m_regangmi_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_regangmi(:,is) = des_m_regangmi(:,is) + H_s2m_regangmi;

        %EIM
        H_s2m_eim = amp_src(is).*exp(1i*prm.k(kk)*(m_eim_pos{kk}(:,1)*cos(theta_s(is))+m_eim_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_eim(:,is) = des_m_eim(:,is) + H_s2m_eim;
    end

    %% Pressure matching

    %Reg
    drv_reg = pinv(H_l2m_reg)*des_m_reg;
    cond_db_reg(kk) = 20*log10(cond(H_l2m_reg));

    %Reg+EqAng
    drv_regang = pinv(H_l2m_regang)*des_m_regang;
    cond_db_regang(kk) = 20*log10(cond(H_l2m_regang));

    %Det
    drv_det = pinv(H_l2m_det)*des_m_det;
    cond_db_det(kk) = 20*log10(cond(H_l2m_det));

    %MI
    drv_mi = pinv(H_l2m_mi)*des_m_mi;
    cond_db_mi(kk) = 20*log10(cond(H_l2m_mi));

    %Det+EqAng 
    drv_regangdet = pinv(H_l2m_regangdet)*des_m_regangdet;
    cond_db_regangdet(kk) = 20*log10(cond(H_l2m_regangdet));

    %MI+EqAng
    drv_regangmi = pinv(H_l2m_regangmi)*des_m_regangmi;
    cond_db_regangmi(kk) = 20*log10(cond(H_l2m_regangmi));

    %EIM
    drv_eim = pinv(H_l2m_eim)*des_m_eim;
    cond_db_eim(kk) = 20*log10(cond(H_l2m_eim));

    %% Reproduction

    fprintf('Reproduction...\n');

    amp_src_mat = ones(rep_prm.num,1)*amp_src;
    ideal = amp_src_mat.*exp(1i*prm.k(kk)*(rep_prm.pos(:,1)*cos(theta_s)+rep_prm.pos(:,2)*sin(theta_s)));

    H_l2r_reg = H_lc2r(:,l_reg_idx{kk},kk);
    H_l2r_regang = H_lc2r(:,l_regang_idx{kk},kk);
    H_l2r_det = H_lc2r(:,l_reg_idx{kk},kk);
    H_l2r_mi = H_lc2r(:,l_reg_idx{kk},kk);
    H_l2r_regangdet = H_lc2r(:,l_regang_idx{kk},kk);
    H_l2r_regangmi = H_lc2r(:,l_regang_idx{kk},kk);
    H_l2r_eim = H_lc2r(:,l_eim_idx{kk},kk);

    %Synthesized pressure distribution
    syn_reg = H_l2r_reg*drv_reg;
    syn_regang = H_l2r_regang*drv_regang;
    syn_det = H_l2r_det*drv_det;
    syn_mi = H_l2r_mi*drv_mi;
    syn_regangdet = H_l2r_regangdet*drv_regangdet;
    syn_regangmi = H_l2r_regangmi*drv_regangmi;
    syn_eim = H_l2r_eim*drv_eim;

    %% Evaluation

    %SDR
    sdrr_reg_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_reg(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_regang_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_regang(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_det_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_det(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_mi_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_mi(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_regangdet_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_regangdet(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_regangmi_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_regangmi(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_eim_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_eim(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));

    sdrr_reg(kk) = mean(sdrr_reg_mat);
    sdrr_regang(kk) = mean(sdrr_regang_mat);
    sdrr_det(kk) = mean(sdrr_det_mat);
    sdrr_mi(kk) = mean(sdrr_mi_mat);
    sdrr_regangdet(kk) = mean(sdrr_regangdet_mat);
    sdrr_regangmi(kk) = mean(sdrr_regangmi_mat);
    sdrr_eim(kk) = mean(sdrr_eim_mat);

    fprintf('[SDR] Reg: %f, Reg+EqAng: %f, Det: %f, MI: %f, Det+EqAng: %f, MI+EqAng: %f, EIM: %f\n',sdrr_reg(kk),sdrr_regang(kk),sdrr_det(kk),sdrr_mi(kk),sdrr_regangdet(kk),sdrr_regangmi(kk),sdrr_eim(kk));
    fprintf('[Cond num] Reg: %f, Reg+EqAng: %f, Det: %f, MI: %f, Det+EqAng: %f, MI+EqAng: %f, EIM: %f\n',cond_db_reg(kk),cond_db_regang(kk),cond_db_det(kk),cond_db_mi(kk),cond_db_regangdet(kk),cond_db_regangmi(kk),cond_db_eim(kk));

    %% Draw figures

    if ismember(prm.freq(kk), figsave_freq) == 1
        idx_plt = 220; %Index for plot

        ideal_plt = reshape(ideal(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_reg_plt = reshape(syn_reg(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_regang_plt = reshape(syn_regang(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_det_plt = reshape(syn_det(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_mi_plt = reshape(syn_mi(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_regangdet_plt = reshape(syn_regangdet(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_regangmi_plt = reshape(syn_regangmi(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_eim_plt = reshape(syn_eim(:,idx_plt),rep_prm.Nx,rep_prm.Ny);

        err_reg_plt = 10*log10((abs(syn_reg_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_regang_plt = 10*log10((abs(syn_regang_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_det_plt = 10*log10((abs(syn_det_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_mi_plt = 10*log10((abs(syn_mi_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_regangdet_plt = 10*log10((abs(syn_regangdet_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_regangmi_plt = 10*log10((abs(syn_regangmi_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_eim_plt = 10*log10((abs(syn_eim_plt-ideal_plt).^2)./(abs(ideal_plt).^2));

        sdrr_reg_plt = sdrr_reg_mat(idx_plt);
        sdrr_regang_plt = sdrr_regang_mat(idx_plt);
        sdrr_det_plt = sdrr_det_mat(idx_plt);
        sdrr_mi_plt = sdrr_mi_mat(idx_plt);
        sdrr_regangdet_plt = sdrr_regangdet_mat(idx_plt);
        sdrr_regangmi_plt = sdrr_regangmi_mat(idx_plt);
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

        %Reg+EqAng
        fig(13)=figure(13);
        set(fig(13),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_regang_plt).');
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

        %Det+EqAng
        fig(14)=figure(14);
        set(fig(14),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_regangdet_plt).');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %MI+EqAng
        fig(17)=figure(17);
        set(fig(17),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_regangmi_plt).');
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

        %Reg+EqAng
        fig(22)=figure(22);
        set(fig(22),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_regang_plt.');
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

        %Det+EqAng
        fig(23)=figure(23);
        set(fig(23),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_regangdet_plt.');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %MI+EqAng
        fig(26)=figure(26);
        set(fig(26),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_regangmi_plt.');
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

        drawnow;
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
plot(prm.freq,sdrr_det,'-^','Color',purple,'MarkerSize',6);
plot(prm.freq,sdrr_regangdet,'-.*','Color',yellow,'MarkerSize',6);
plot(prm.freq,sdrr_eim,'-+','Color',brown,'MarkerSize',6);
hold off;
grid on;
xlim([20,1600]);
ylim([0,70]);
xlabel('Frequency (Hz)'); ylabel('SDR (dB)');
legend('Reg','Det','Det+EqAng','EIM','Location','southwest');

fig(3)=figure(3);
set(fig(3),'Position',fig_pos);
hold on;
plot(prm.freq,sdrr_reg,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,sdrr_mi,'--v','Color',green,'MarkerSize',6);
plot(prm.freq,sdrr_regangmi,'-.s','Color',sky,'MarkerSize',6);
plot(prm.freq,sdrr_eim,'-+','Color',brown,'MarkerSize',6);
hold off;
grid on;
xlim([20,1600]);
ylim([0,70]);
xlabel('Frequency (Hz)'); ylabel('SDR (dB)');
legend('Reg','MI','MI+EqAng','EIM','Location','southwest');

%Condition number
fig(4)=figure(4);
set(fig(4),'Position',fig_pos);
hold on;
plot(prm.freq,cond_db_reg,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,cond_db_det,'-^','Color',purple,'MarkerSize',6);
plot(prm.freq,cond_db_regangdet,'-.*','Color',yellow,'MarkerSize',6);
plot(prm.freq,cond_db_eim,'-+','Color',brown,'MarkerSize',6);
hold off;
grid on;
xlim([20,1600]);
ylim([40,150]);
xlabel('Frequency (Hz)'); ylabel('Condition number (dB)');
legend('Reg','Det','Det+EqAng','EIM','Location','northeast');

fig(5)=figure(5);
set(fig(5),'Position',fig_pos);
hold on;
plot(prm.freq,cond_db_reg,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,cond_db_mi,'--v','Color',green,'MarkerSize',6);
plot(prm.freq,cond_db_regangmi,'-.s','Color',sky,'MarkerSize',6);
plot(prm.freq,cond_db_eim,'-+','Color',brown,'MarkerSize',6);
hold off;
grid on;
xlim([20,1600]);
ylim([40,150]);
xlabel('Frequency (Hz)'); ylabel('Condition number (dB)');
legend('Reg','MI','MI+EqAng','EIM','Location','northeast');

if figsave_flg==1
    fname = cell(2,1);

    fname{1} = sprintf('%s/fig15a',dir_results);
    fname{2} = sprintf('%s/fig15b',dir_results);
    save_figures(fig(2:3), fname, fig_type);
end
