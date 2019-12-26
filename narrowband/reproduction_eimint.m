% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

close all;

rnd_rep_flg = 0; %Flag for evaluating random sampling of internal points

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
sdrr_regchief = zeros(prm.freq_len,1);
sdrr_regrndi = zeros(prm.freq_len,m_rndi_ntrial);
sdrr_regeim = zeros(prm.freq_len,1);
sdrr_regeimin = zeros(prm.freq_len,m_eimi_num_max);
sdrr_eim = zeros(prm.freq_len,1);

%Condition number in dB
cond_db_reg = zeros(prm.freq_len,1);
cond_db_regchief = zeros(prm.freq_len,m_chief_num);
cond_db_regrndi = zeros(prm.freq_len,m_rndi_ntrial);
cond_db_regrndi_best = zeros(prm.freq_len,1);
cond_db_regeim = zeros(prm.freq_len,1);
cond_db_regeimin = zeros(prm.freq_len,m_eimi_num_max);
cond_db_eim = zeros(prm.freq_len,1);

if rnd_rep_flg == 0
    %Load results of random sampling
    if rev_prm.flg == 1
        fname = sprintf('../data/pos_int/abs%.2f_x%d_y%d/int%d_trial%d/result_rep.mat',rev_prm.abs_ratio,des_prm.len_x*100, des_prm.len_y*100,m_int_num,m_rndi_ntrial);
    else
        fname = sprintf('../data/pos_int/freefield_x%d_y%d/int%d_trial%d/result_rep.mat',des_prm.len_x*100, des_prm.len_y*100,m_int_num,m_rndi_ntrial);
    end
    load(fname,'sdrr_regrndi','cond_db_regrndi');
end

for kk=1:prm.freq_len
    fprintf('freq: %d\n',prm.freq(kk));
    
    %% Random sampling

    fname_rndi_pos = sprintf('%s/pos_rndi_f%d.mat',dirname_rndi_pos,prm.freq(kk));
    load(fname_rndi_pos,'m_regrndi_num','m_regrndi_idx','m_regrndi_pos');

    amp_src_mat = ones(rep_prm.num,1)*amp_src;
    ideal = amp_src_mat.*exp(1i*prm.k(kk)*(rep_prm.pos(:,1)*cos(theta_s)+rep_prm.pos(:,2)*sin(theta_s)));
    H_l2r_reg = H_lc2r(:,l_reg_idx{kk},kk);

    if rnd_rep_flg == 1
        sdrr_regrndi_mat = zeros(Nsrc,m_rndi_ntrial);

        for ii=1:m_rndi_ntrial
            H_l2m_regrndi = H_lc2d(m_regrndi_idx(:,ii),l_reg_idx{kk},kk);
            des_m_regrndi = zeros(m_regrndi_num,Nsrc);
            for is=1:Nsrc
                H_s2m_regrndi = amp_src(is).*exp(1i*prm.k(kk)*(m_regrndi_pos(:,1,ii)*cos(theta_s(is))+m_regrndi_pos(:,2,ii)*sin(theta_s(is))));
                des_m_regrndi(:,is) = des_m_regrndi(:,is) + H_s2m_regrndi;
            end
            drv_regrndi = pinv(H_l2m_regrndi)*des_m_regrndi;
            cond_db_regrndi(kk,ii) = 20*log10(cond(H_l2m_regrndi));

            syn_regrndi = H_l2r_reg*drv_regrndi;
            sdrr_regrndi_mat(:,ii) = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_regrndi(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
            sdrr_regrndi(kk,ii) = mean(sdrr_regrndi_mat(:,ii));
        end
    end

    [~,sdrr_regrnd_max_idx]=max(sdrr_regrndi(kk,:));

    %Best result
    H_l2m_regrndi = H_lc2d(m_regrndi_idx(:,sdrr_regrnd_max_idx),l_reg_idx{kk},kk);
    cond_db_regrndi_best(kk) = 20*log10(cond(H_l2m_regrndi));
    des_m_regrndi = zeros(m_regrndi_num,Nsrc);
    for is=1:Nsrc
        H_s2m_regrndi = amp_src(is).*exp(1i*prm.k(kk)*(m_regrndi_pos(:,1,sdrr_regrnd_max_idx)*cos(theta_s(is))+m_regrndi_pos(:,2,sdrr_regrnd_max_idx)*sin(theta_s(is))));
        des_m_regrndi(:,is) = des_m_regrndi(:,is) + H_s2m_regrndi;
    end
    drv_regrndi = pinv(H_l2m_regrndi)*des_m_regrndi;
    syn_regrndi = H_l2r_reg*drv_regrndi;

    %% Calculate transfer functions

    H_l2m_reg = H_lc2d(m_reg_idx{kk},l_reg_idx{kk},kk);
    H_l2m_regchief = H_lc2d(m_regchief_idx{kk},l_reg_idx{kk},kk);
    H_l2m_regeim = H_lc2d(m_regeim_idx{kk},l_reg_idx{kk},kk);
    H_l2m_regeimin = cell(m_eimi_num_max,1);
    for ii=1:m_eimi_num_max
        H_l2m_regeimin{ii} = H_lc2d(m_regeimin_idx{kk,ii},l_reg_idx{kk},kk);
    end
    H_l2m_eim = H_lc2d(m_eim_idx{kk},l_eim_idx{kk},kk);

    %% Desired sound field 

    des = zeros(des_prm.num,Nsrc);
    des_m_reg = zeros(m_reg_num(kk),Nsrc);
    des_m_regchief = zeros(m_regchief_num(kk),Nsrc);
    des_m_regeim = zeros(m_regeim_num(kk),Nsrc);
    des_m_regeimin = cell(m_eimi_num_max,1);
    for ii=1:m_eimi_num_max
        des_m_regeimin{ii} = zeros(m_regeimin_num(kk,ii),Nsrc);
    end
    des_m_eim = zeros(m_eim_num(kk),Nsrc);

    for is=1:Nsrc
        %Target
        H_s2d = amp_src(is).*exp(1i*prm.k(kk)*(des_prm.pos(:,1)*cos(theta_s(is))+des_prm.pos(:,2)*sin(theta_s(is))));
        des(:,is) = des(:,is) + H_s2d;

        %Reg
        H_s2m_reg = amp_src(is).*exp(1i*prm.k(kk)*(m_reg_pos{kk}(:,1)*cos(theta_s(is))+m_reg_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_reg(:,is) = des_m_reg(:,is) + H_s2m_reg;

        %Reg+CHIEF
        H_s2m_regchief = amp_src(is).*exp(1i*prm.k(kk)*(m_regchief_pos{kk}(:,1)*cos(theta_s(is))+m_regchief_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_regchief(:,is) = des_m_regchief(:,is) + H_s2m_regchief;

        %Reg+EIM
        H_s2m_regeim = amp_src(is).*exp(1i*prm.k(kk)*(m_regeim_pos{kk}(:,1)*cos(theta_s(is))+m_regeim_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_regeim(:,is) = des_m_regeim(:,is) + H_s2m_regeim;

        %Reg+EIMi
        for ii=1:m_eimi_num_max
            H_s2m_regeimin = amp_src(is).*exp(1i*prm.k(kk)*(m_regeimin_pos{kk,ii}(:,1)*cos(theta_s(is))+m_regeimin_pos{kk,ii}(:,2)*sin(theta_s(is))));
            des_m_regeimin{ii}(:,is) = des_m_regeimin{ii}(:,is) + H_s2m_regeimin;
        end

        %EIM
        H_s2m_eim = amp_src(is).*exp(1i*prm.k(kk)*(m_eim_pos{kk}(:,1)*cos(theta_s(is))+m_eim_pos{kk}(:,2)*sin(theta_s(is))));
        des_m_eim(:,is) = des_m_eim(:,is) + H_s2m_eim;
    end

    %% Pressure matching

    %Reg
    drv_reg = pinv(H_l2m_reg)*des_m_reg;
    cond_db_reg(kk) = 20*log10(cond(H_l2m_reg));

    %Reg+CHIEF
    drv_regchief = pinv(H_l2m_regchief)*des_m_regchief;
    cond_db_regchief(kk) = 20*log10(cond(H_l2m_regchief));

    %Reg+EIM
    drv_regeim = pinv(H_l2m_regeim)*des_m_regeim;
    cond_db_regeim(kk) = 20*log10(cond(H_l2m_regeim));

    %Reg+EIMi
    drv_regeimin = cell(m_eimi_num_max,1);
    for ii=1:m_eimi_num_max
        drv_regeimin{ii} = pinv(H_l2m_regeimin{ii})*des_m_regeimin{ii};
        cond_db_regeimin(kk,ii) = 20*log10(cond(H_l2m_regeimin{ii}));
    end

    %EIM
    drv_eim = pinv(H_l2m_eim)*des_m_eim;
    cond_db_eim(kk) = 20*log10(cond(H_l2m_eim));

    %% Reproduction

    fprintf('Reproduction...\n');

    H_l2r_reg = H_lc2r(:,l_reg_idx{kk},kk);
    H_l2r_eim = H_lc2r(:,l_eim_idx{kk},kk);

    %Synthesized pressure distribution
    syn_reg = H_l2r_reg*drv_reg;

    syn_regchief = H_l2r_reg*drv_regchief;

    syn_regeim = H_l2r_reg*drv_regeim;

    syn_regeimin = cell(m_eimi_num_max,1);
    for ii=1:m_eimi_num_max
        syn_regeimin{ii} = H_l2r_reg*drv_regeimin{ii};
    end

    syn_eim = H_l2r_eim*drv_eim;

    %% Evaluation

    %SDR
    sdrr_reg_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_reg(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));

    sdrr_regchief_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_regchief(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    
    sdrr_regeim_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_regeim(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    sdrr_regeimin_mat = zeros(Nsrc,m_eimi_num_max);
    for ii=1:m_eimi_num_max
        sdrr_regeimin_mat(:,ii) = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_regeimin{ii}(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));
    end

    sdrr_eim_mat = 10*log10(sum(abs(ideal(rep_prm.idx_des,:)).^2,1)./sum(abs(syn_eim(rep_prm.idx_des,:)-ideal(rep_prm.idx_des,:)).^2,1));

    sdrr_reg(kk) = mean(sdrr_reg_mat);
    sdrr_regchief(kk) = mean(sdrr_regchief_mat);

    sdrr_regeim(kk) = mean(sdrr_regeim_mat);
    for ii=1:m_eimi_num_max
        sdrr_regeimin(kk,ii) = mean(sdrr_regeimin_mat(:,ii));
    end
    sdrr_eim(kk) = mean(sdrr_eim_mat);

    fprintf('[SDR] Reg: %f, Reg+CHIEF: %f, Random(best): %f, Reg+EIM: %f, Reg+EIMi: %f, EIM: %f\n',sdrr_reg(kk),sdrr_regchief(kk),sdrr_regrndi(kk,sdrr_regrnd_max_idx),sdrr_regeim(kk),sdrr_regeimin(kk,m_int_num),sdrr_eim(kk));
    fprintf('[Cond num] Reg: %f, Reg+CHIEF: %f, Random(best): %f, Reg+EIM: %f, Reg+EIMi: %f, EIM: %f\n',cond_db_reg(kk),cond_db_regchief(kk),cond_db_regrndi(kk,sdrr_regrnd_max_idx),cond_db_regeim(kk),cond_db_regeimin(kk,m_int_num),cond_db_eim(kk));

    %% Draw figures

    if ismember(prm.freq(kk), figsave_freq) == 1
        idx_plt = 220; %Index for plot

        ideal_plt = reshape(ideal(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_reg_plt = reshape(syn_reg(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_regchief_plt = reshape(syn_regchief(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_regrndi_plt = reshape(syn_regrndi(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_regeim_plt = reshape(syn_regeim(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_regeimin_plt = reshape(syn_regeimin{m_int_num}(:,idx_plt),rep_prm.Nx,rep_prm.Ny);
        syn_eim_plt = reshape(syn_eim(:,idx_plt),rep_prm.Nx,rep_prm.Ny);

        err_reg_plt = 10*log10((abs(syn_reg_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_regchief_plt = 10*log10((abs(syn_regchief_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_regrndi_plt = 10*log10((abs(syn_regrndi_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_regeim_plt = 10*log10((abs(syn_regeim_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_regeimin_plt = 10*log10((abs(syn_regeimin_plt-ideal_plt).^2)./(abs(ideal_plt).^2));
        err_eim_plt = 10*log10((abs(syn_eim_plt-ideal_plt).^2)./(abs(ideal_plt).^2));

        sdrr_reg_plt = sdrr_reg_mat(idx_plt);
        sdrr_regchief_plt = sdrr_regchief_mat(idx_plt);
        sdrr_regrndi_plt = 10*log10(sum(abs(ideal(rep_prm.idx_des,idx_plt)).^2,1)./sum(abs(syn_regrndi(rep_prm.idx_des,idx_plt)-ideal(rep_prm.idx_des,idx_plt)).^2,1));
        sdrr_regeim_plt = sdrr_regeim_mat(idx_plt);
        sdrr_regeimin_plt = sdrr_regeimin_mat(idx_plt,m_int_num);
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

        %Reg+CHIEF
        fig(14)=figure(14);
        set(fig(14),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_regchief_plt).');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %Reg+EIMi
        fig(16)=figure(16);
        set(fig(16),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_regeimin_plt).');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %Reg+EIM
        fig(17)=figure(17);
        set(fig(17),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_regeim_plt).');
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

        %Best result of random sampling
        fig(19)=figure(19);
        set(fig(19),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,real(syn_regrndi_plt).');
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
        fig(31)=figure(31);
        set(fig(31),'Position',fig_pos);
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

        %Reg+CHIEF
        fig(33)=figure(33);
        set(fig(33),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_regchief_plt.');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %Reg+EIMi
        fig(35)=figure(35);
        set(fig(35),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_regeimin_plt.');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        %Reg+EIM
        fig(36)=figure(36);
        set(fig(36),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_regeim_plt.');
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
        fig(37)=figure(37);
        set(fig(37),'Position',fig_pos);
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

        %Best result of random sampling
        fig(38)=figure(38);
        set(fig(38),'Position',fig_pos);
        hold on;
        imagesc(xx,yy,err_regrndi_plt.');
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        hold off;
        axis tight;
        axis equal;
        caxis(z_range);
        colormap(flipud(pink));
        colorbar;
        set(gca,'Ydir','normal');
        xlabel('x (m)'); ylabel('y (m)');

        % Best placement of random sampling
        fig_pos = [0, 0, 300, 350];

        fig(51)=figure(51);
        set(fig(51),'Position',fig_pos);
        hold on;
        plot(m_regrndi_pos(:,1,sdrr_regrnd_max_idx),m_regrndi_pos(:,2,sdrr_regrnd_max_idx),'xk','MarkerSize',4); %Microphones
        plot(lc_prm.pos(:,1),lc_prm.pos(:,2),'.k','MarkerSize',1);
        plot(l_reg_pos{kk}(:,1),l_reg_pos{kk}(:,2),'ok','MarkerSize',4,'MarkerFaceColor','k'); %Loudspeakers
        rectangle('Position',[-des_prm.len_x/2, -des_prm.len_y/2, des_prm.len_x, des_prm.len_y],'LineWidth',1,'LineStyle','--');
        if rev_prm.flg == 1; plot([rev_prm.pos_bl(1), rev_prm.pos_br(1), rev_prm.pos_tr(1), rev_prm.pos_tl(1), rev_prm.pos_bl(1)], [rev_prm.pos_bl(2), rev_prm.pos_br(2), rev_prm.pos_tr(2), rev_prm.pos_tl(2), rev_prm.pos_bl(2)],'-k','LineWidth',2); end;
        hold off;
        xlabel('x (m)'); ylabel('y (m)');
        axis equal;
        if rev_prm.flg == 1
            xlim([rev_prm.pos_bl(1)*1.1,rev_prm.pos_br(1)*1.1]);
            ylim([rev_prm.pos_bl(2)*1.1,rev_prm.pos_tl(2)*1.1]);
        else 
            xlim([min(lc_prm.pos(:,1))*1.1,max(lc_prm.pos(:,1))*1.1]);
            ylim([min(lc_prm.pos(:,2))*1.1,max(lc_prm.pos(:,2))*1.1]);
        end

        drawnow;
        
        if figsave_flg==1
            fname = cell(1,1);

            fname{1} = sprintf('%s/fig11c',dir_results);
            save_figures(fig(51), fname, fig_type);
        end
    end

    clear 'm_regrndi_num' 'm_regrndi_idx' 'm_regrndi_pos';
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

sdrr_regrndi_tmean = mean(sdrr_regrndi,2);
sdrr_regrndi_tmax = max(sdrr_regrndi,[],2);
sdrr_regrndi_tmin = min(sdrr_regrndi,[],2);

%SDR
fig(2)=figure(2);
set(fig(2),'Position',fig_pos);
hold on;
plot(prm.freq,sdrr_reg,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,sdrr_regchief,'-.*','Color',yellow,'MarkerSize',6);
plot(prm.freq,sdrr_regeimin(:,6),'--v','Color',green,'MarkerSize',6);
plot(prm.freq,sdrr_eim,'-+','Color',brown,'MarkerSize',6);
h = area(prm.freq,[sdrr_regrndi_tmin,sdrr_regrndi_tmax-sdrr_regrndi_tmin]);
set(h(1),'FaceColor','None');
set(h(2),'FaceColor',[0.0,0.1,0.1],'FaceAlpha',0.1);
set(h(1),'EdgeColor','None');
set(h(2),'EdgeColor','None');
hold off;
grid on;
xlim([20,1600]);
ylim([0,60]);
xlabel('Frequency (Hz)'); ylabel('SDR (dB)');
legend('Reg','Reg+CHIEF','Reg+EIMi','EIM','Location','southwest');

fig(5)=figure(5);
set(fig(5),'Position',fig_pos);
hold on;
plot(prm.freq,sdrr_reg,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,sdrr_regeimin(:,2),'--^','Color',red,'MarkerSize',6);
plot(prm.freq,sdrr_regeimin(:,6),'--v','Color',green,'MarkerSize',6);
plot(prm.freq,sdrr_regeimin(:,10),'-.*','Color',yellow,'MarkerSize',6);
plot(prm.freq,sdrr_eim(:,1),'-+','Color',brown,'MarkerSize',6);
h = area(prm.freq,[sdrr_regrndi_tmin,sdrr_regrndi_tmax-sdrr_regrndi_tmin]);
set(h(1),'FaceColor','None');
set(h(2),'FaceColor',[0.0,0.1,0.1],'FaceAlpha',0.1);
set(h(1),'EdgeColor','None');
set(h(2),'EdgeColor','None');
hold off;
grid on;
xlim([20,1600]);
ylim([0,60]);
xlabel('Frequency (Hz)'); ylabel('SDR (dB)');
legend('Reg','Reg+EIMi(2)','Reg+EIMi(6)','Reg+EIMi(10)','EIM','Location','southwest');

fig(6)=figure(6);
set(fig(6),'Position',fig_pos);
hold on;
plot(prm.freq,sdrr_reg,'-o','Color',blue,'MarkerSize',6);
plot(prm.freq,sdrr_regeim,'-.x','Color',red,'MarkerSize',6);
plot(prm.freq,sdrr_regeimin(:,m_int_num),'--v','Color',green,'MarkerSize',6);
plot(prm.freq,sdrr_eim,'-+','Color',brown,'MarkerSize',6);
h = area(prm.freq,[sdrr_regrndi_tmin,sdrr_regrndi_tmax-sdrr_regrndi_tmin]);
set(h(1),'FaceColor','None');
set(h(2),'FaceColor',[0.0,0.1,0.1],'FaceAlpha',0.1);
set(h(1),'EdgeColor','None');
set(h(2),'EdgeColor','None');
hold off;
grid on;
xlim([20,1600]);
ylim([0,60]);
xlabel('Frequency (Hz)'); ylabel('SDR (dB)');
legend('Reg','Reg+EIM','Reg+EIMi','EIM','Location','southwest');

if figsave_flg==1
    fname = cell(3,1);

    fname{1} = sprintf('%s/fig10',dir_results);
    fname{2} = sprintf('%s/fig12',dir_results);
    fname{3} = sprintf('%s/fig13',dir_results);
    save_figures([fig(2),fig(5:6)], fname, fig_type);
end
