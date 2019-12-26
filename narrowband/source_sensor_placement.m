% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

clear variables;
close all;

rng(0);

addpath('../routine');
addpath('../routine_external');

set(0,'defaultAxesFontSize',18);
set(0,'defaultAxesFontName','Times');
set(0,'defaultTextFontSize',18);
set(0,'defaultTextFontName','Times');

ColorList;

figsave_freq = 800; %Frequency of generating figures
figsave_flg = 1; %Flag for saving figures
dir_results = '../results'; %Directory to save figures
fig_type = 'png'; %Figure type (png/pdf/epsc/fig)
pos_rnd_gen_flg = 0; %Flag for regenerating Rand

%% Parameters

%Sound velocity (m/s)
prm.c=340.0;

%Frequency (Hz)
prm.freq = 20:20:1600;
prm.freq_len = length(prm.freq);

%Wave number (1/m)
prm.k=2*pi*prm.freq./prm.c;

%% Parameters: reverberation

rev_prm.flg = 1; %0: free field, 1: reverberant

if rev_prm.flg == 1
    rev_prm.abs_ratio = 0.5; %Absorption ratio of walls
    fprintf('Absorption ratio: %f\n',rev_prm.abs_ratio);

    %Acoustic impedance ratio
    rev_prm.imp_ratio = (2-rev_prm.abs_ratio+2*sqrt(1-rev_prm.abs_ratio))/rev_prm.abs_ratio;

    %Room geometry
    rev_prm.pos_bl = [-2.5, -3.0]; %bottom left
    rev_prm.pos_br = [2.5, -3.0]; %bottom right
    rev_prm.pos_tr = [2.5, 3.0]; %top right
    rev_prm.pos_tl = [-2.5, 4.0]; %top left

    rev_prm.len_l = sqrt((rev_prm.pos_tl(1)-rev_prm.pos_bl(1))^2+(rev_prm.pos_tl(2)-rev_prm.pos_bl(2))^2); %left
    rev_prm.len_r = sqrt((rev_prm.pos_tr(1)-rev_prm.pos_br(1))^2+(rev_prm.pos_tr(2)-rev_prm.pos_br(2))^2); %right
    rev_prm.len_t = sqrt((rev_prm.pos_tl(1)-rev_prm.pos_tr(1))^2+(rev_prm.pos_tl(2)-rev_prm.pos_tr(2))^2); %top
    rev_prm.len_b = sqrt((rev_prm.pos_bl(1)-rev_prm.pos_br(1))^2+(rev_prm.pos_bl(2)-rev_prm.pos_br(2))^2); %bottom
end

%% Parameters: Reproduced area

%Size of visualizing rectangular region
rep_prm.len_x = 1.5;
rep_prm.len_y = 1.5;

%Interval
rep_prm.dx = 0.02;
rep_prm.dy = 0.02;

%Number of samples
rep_prm.Nx=round(rep_prm.len_x/rep_prm.dx);
rep_prm.Ny=round(rep_prm.len_y/rep_prm.dy);

%Position
rep_prm.x = (((0:rep_prm.Nx-1)-rep_prm.Nx/2).*rep_prm.dx)';
rep_prm.y = (((0:rep_prm.Ny-1)-rep_prm.Ny/2).*rep_prm.dy)';

x_r_vec = reshape(rep_prm.x*ones(1,rep_prm.Ny),rep_prm.Nx*rep_prm.Ny,1);
y_r_vec = reshape(ones(rep_prm.Nx,1)*rep_prm.y.',rep_prm.Nx*rep_prm.Ny,1);

rep_prm.pos = [x_r_vec, y_r_vec];
rep_prm.num = size(rep_prm.pos,1);

%% Parameters: Desired area

%Size of rectangular target region
des_prm.len_x = 0.8;
des_prm.len_y = 1.0;

%Interval
des_prm.dx = 0.04;
des_prm.dy = 0.04;

%Number of samples
des_prm.Nx=round(des_prm.len_x/des_prm.dx)+1;
des_prm.Ny=round(des_prm.len_y/des_prm.dy)+1;

%Position
des_prm.x = (((0:des_prm.Nx-1)-des_prm.Nx/2).*des_prm.dx+des_prm.dx/2)';
des_prm.y = (((0:des_prm.Ny-1)-des_prm.Ny/2).*des_prm.dy+des_prm.dy/2)';

x_d_vec = reshape(des_prm.x*ones(1,des_prm.Ny),des_prm.Nx*des_prm.Ny,1);
y_d_vec = reshape(ones(des_prm.Nx,1)*des_prm.y.',des_prm.Nx*des_prm.Ny,1);

des_prm.pos = [x_d_vec, y_d_vec];
des_prm.x = des_prm.pos(:,1);
des_prm.y = des_prm.pos(:,2);
des_prm.num = size(des_prm.pos,1);

%Indexes of target region inside reproduced area
rep_prm.idx_des = rep_prm.pos(:,1)<=des_prm.len_x/2 & rep_prm.pos(:,1)>=-des_prm.len_x/2 & rep_prm.pos(:,2)<=des_prm.len_y/2 & rep_prm.pos(:,2)>=-des_prm.len_y/2;

%% Parameters: Loudspeaker candidates

%Size of rectangular loudspeaker candidates
lc_prm.len_x = 2.4;
lc_prm.len_y = 2.8;

%Position of center
lc_prm.shift_x = -0.1;
lc_prm.shift_y = -0.2;

%Number of candidate points
lc_prm.num = 256;

%Positions on boundary of rectangular region
lc_prm.pos = rect_perim(lc_prm.num, lc_prm.len_x, lc_prm.len_y, lc_prm.shift_x, lc_prm.shift_y);

%% Calculate transfer functions of loudspeakers

fprintf('Calculating transfer functions...\n');

pos_rd = [rep_prm.pos; des_prm.pos];
pos_rd_num = size(pos_rd,1);

H_lc2rd = zeros(pos_rd_num, lc_prm.num, prm.freq_len);
if rev_prm.flg==0
    %Freefield
    for kk=1:prm.freq_len
        H_lc2rd(:,:,kk) = green2d(repmat(pos_rd(:,1),1,lc_prm.num),repmat(pos_rd(:,2),1,lc_prm.num),repmat(lc_prm.pos(:,1).',pos_rd_num,1),repmat(lc_prm.pos(:,2).',pos_rd_num,1),prm.k(kk));
    end
else
    %Reverberant (load files)
    for kk=1:prm.freq_len
        fname_imp_lc2rd_f = sprintf('../data/imp_lc2rd_abs%.2f_l%d_%dx%d_r%d_%dx%d_d%d_%dx%d/imp_lc2rd_f%d.mat',rev_prm.abs_ratio,lc_prm.num,lc_prm.len_x*100,lc_prm.len_y*100,rep_prm.num,rep_prm.len_x*100,rep_prm.len_y*100,des_prm.num,des_prm.len_x*100,des_prm.len_y*100,prm.freq(kk));
        load(fname_imp_lc2rd_f,'H_lc2rd_f');
        H_lc2rd(:,:,kk) = H_lc2rd_f;
    end
end

%Transfer functions from loudspeaker candidates to reproduced area
H_lc2r = H_lc2rd(1:rep_prm.num,:,:);
%Transfer functions from loudspeaker candidates to desired area
H_lc2d = H_lc2rd(rep_prm.num+1:rep_prm.num+des_prm.num,:,:);

clear H_lc2rd;

%% Memory allocation

%Reg
l_reg_num = zeros(prm.freq_len,1);
l_reg_idx = cell(prm.freq_len,1);
l_reg_pos = cell(prm.freq_len,1);
m_reg_num = zeros(prm.freq_len,1);
m_reg_idx = cell(prm.freq_len,1);
m_reg_pos = cell(prm.freq_len,1);

%Rand
fname_rnd = sprintf('./pos_rnd_x%d_y%d.mat',des_prm.len_x*100, des_prm.len_y*100);

if pos_rnd_gen_flg == 1
    l_rnd_num = zeros(prm.freq_len,1);
    l_rnd_idx = cell(prm.freq_len,1);
    l_rnd_pos = cell(prm.freq_len,1);
    m_rnd_num = zeros(prm.freq_len,1);
    m_rnd_idx = cell(prm.freq_len,1);
    m_rnd_pos = cell(prm.freq_len,1);
else
    load(fname_rnd);
end

%GSO
l_gso_num = zeros(prm.freq_len,1);
l_gso_idx = cell(prm.freq_len,1);
l_gso_pos = cell(prm.freq_len,1);
m_gso_num = zeros(prm.freq_len,1);
m_gso_idx = cell(prm.freq_len,1);
m_gso_pos = cell(prm.freq_len,1);

%Det
l_det_num = zeros(prm.freq_len,1);
l_det_idx = cell(prm.freq_len,1);
l_det_pos = cell(prm.freq_len,1);
m_det_num = zeros(prm.freq_len,1);
m_det_idx = cell(prm.freq_len,1);
m_det_pos = cell(prm.freq_len,1);

%MI
l_mi_num = zeros(prm.freq_len,1);
l_mi_idx = cell(prm.freq_len,1);
l_mi_pos = cell(prm.freq_len,1);
m_mi_num = zeros(prm.freq_len,1);
m_mi_idx = cell(prm.freq_len,1);
m_mi_pos = cell(prm.freq_len,1);

%FS
l_fs_num = zeros(prm.freq_len,1);
l_fs_idx = cell(prm.freq_len,1);
l_fs_pos = cell(prm.freq_len,1);
m_fs_num = zeros(prm.freq_len,1);
m_fs_idx = cell(prm.freq_len,1);
m_fs_pos = cell(prm.freq_len,1);

%EIM
l_eim_num = zeros(prm.freq_len,1);
l_eim_idx = cell(prm.freq_len,1);
l_eim_pos = cell(prm.freq_len,1);
m_eim_num = zeros(prm.freq_len,1);
m_eim_idx = cell(prm.freq_len,1);
m_eim_pos = cell(prm.freq_len,1);

%% Source/sensor placement (loudspeakers - microphones)

for kk=1:prm.freq_len
    fprintf('freq: %d\n',prm.freq(kk));

    %% Empirical interpolation method (EIM)

    lmp_eim_thr = 1e-2;
    [l_eim_idx{kk}, m_eim_idx{kk}] = lmp_s_eim_th(H_lc2d(:,:,kk), lmp_eim_thr);

    %Loudspeakers
    l_eim_num(kk) = length(l_eim_idx{kk});
    l_eim_pos{kk} = lc_prm.pos(l_eim_idx{kk},:);

    %Microphones
    m_eim_num(kk) = length(m_eim_idx{kk});
    m_eim_pos{kk} = des_prm.pos(m_eim_idx{kk},:);

    %% Regular - Regular (Reg)

    %Loudspeakers
    l_reg_num(kk) = l_eim_num(kk);

    l_reg_pos_tmp = rect_perim(l_reg_num(kk), lc_prm.len_x, lc_prm.len_y, lc_prm.shift_x, lc_prm.shift_y);

    l_reg_idx{kk} = zeros(l_reg_num(kk),1);
    for ii=1:l_reg_num(kk)
        dd = sqrt((l_reg_pos_tmp(ii,1)-lc_prm.pos(:,1)).^2 + (l_reg_pos_tmp(ii,2)-lc_prm.pos(:,2)).^2);
        [~, l_reg_idx{kk}(ii)] = min(dd);
    end

    l_reg_pos{kk} = lc_prm.pos(l_reg_idx{kk},:);

    %Microphones
    m_reg_num(kk) = m_eim_num(kk);

    m_reg_pos_tmp = rect_perim(m_reg_num(kk), des_prm.len_x, des_prm.len_y, 0, 0);

    m_reg_idx{kk} = zeros(m_reg_num(kk),1);
    for ii=1:m_reg_num(kk)
        dd = sqrt((m_reg_pos_tmp(ii,1)-des_prm.pos(:,1)).^2 + (m_reg_pos_tmp(ii,2)-des_prm.pos(:,2)).^2);
        [~, m_reg_idx{kk}(ii)] = min(dd);
    end

    m_reg_pos{kk} = des_prm.pos(m_reg_idx{kk},:);

    %% Random - Random (Rand)

    if pos_rnd_gen_flg == 1
        %Loudspeakers
        l_rnd_num(kk) = l_eim_num(kk);
        l_rnd_idx{kk} = randsample(lc_prm.num,l_rnd_num(kk));

        l_rnd_pos{kk} = lc_prm.pos(l_rnd_idx{kk},:);

        %Microphones
        m_rnd_num(kk) = m_eim_num(kk);
        m_rnd_idx{kk} = randsample(des_prm.num,m_rnd_num(kk));

        m_rnd_pos{kk} = des_prm.pos(m_rnd_idx{kk},:);
    end

    %% Gram-Schmidt Orthogonalization - two layer (GSO)

    %Microphones
    m_gso_num(kk) = m_eim_num(kk);

    num_ext = ceil(m_gso_num(kk)/2);
    num_int = m_gso_num(kk)-num_ext;

    prm.d_2L = 0.04;

    pos_ext = rect_perim(num_ext, des_prm.len_x, des_prm.len_y, 0, 0);
    pos_int = rect_perim(num_int, des_prm.len_x-prm.d_2L*2, des_prm.len_y-prm.d_2L*2, 0, 0);

    m_gso_pos_tmp = [pos_ext; pos_int];

    m_gso_idx{kk} = zeros(m_gso_num(kk),1);
    for ii=1:m_gso_num(kk)
        dd = sqrt((m_gso_pos_tmp(ii,1)-des_prm.pos(:,1)).^2 + (m_gso_pos_tmp(ii,2)-des_prm.pos(:,2)).^2);
        [~, m_gso_idx{kk}(ii)] = min(dd);
    end

    m_gso_pos{kk} = des_prm.pos(m_gso_idx{kk},:);

    %Loudspeakers
    l_gso_num(kk) = l_eim_num(kk);

    %Desired field for initialization
    theta_s = 0;
    des = exp(1i*prm.k(kk)*(des_prm.pos(m_gso_idx{kk},1)*cos(theta_s)+des_prm.pos(m_gso_idx{kk},2)*sin(theta_s)));

    l_gso_idx{kk} = lp_s_gso(H_lc2d(m_gso_idx{kk},:,kk), l_gso_num(kk), des);
    l_gso_pos{kk} = lc_prm.pos(l_gso_idx{kk},:);

    %% Regular - Determinant (Det)

    %Loudspeakers
    l_det_num(kk) = l_reg_num(kk);
    l_det_idx{kk} = l_reg_idx{kk};
    l_det_pos{kk} = l_reg_pos{kk};

    %Microphones
    m_det_num(kk) = m_eim_num(kk);
    m_det_loc_thr = 0.4;
    [~, ~, zast, ~] = mp_s_det_app(H_lc2d(:,l_det_idx{kk},kk), m_det_num(kk));
    [z_loc, ~] =  mp_s_det_locr(H_lc2d(:,l_det_idx{kk},kk), m_det_num(kk), zast, m_det_loc_thr);

    m_det_idx{kk} = find(z_loc);
    m_det_pos{kk} = des_prm.pos(m_det_idx{kk},:);

    %% Regular - Mutual Information (MI)

    %Loudspeakers
    l_mi_num(kk) = l_reg_num(kk);
    l_mi_idx{kk} = l_reg_idx{kk};
    l_mi_pos{kk} = l_reg_pos{kk};

    %Microphones
    m_mi_num(kk) = m_eim_num(kk);
    m_mi_idx{kk} = mp_s_mi(H_lc2d(:,l_mi_idx{kk},kk)*H_lc2d(:,l_mi_idx{kk},kk)',m_mi_num(kk));
    m_mi_pos{kk} = des_prm.pos(m_mi_idx{kk},:);

    %% FrameSense (FS)

    %Loudspeakers
    l_fs_num(kk) = l_eim_num(kk);
    l_fs_idx{kk} = lp_s_fp(H_lc2d(:,:,kk),l_fs_num(kk));
    l_fs_pos{kk} = lc_prm.pos(l_fs_idx{kk},:);

    %Microphones
    m_fs_num(kk) = m_eim_num(kk);
    m_fs_idx{kk} = mp_s_fp(H_lc2d(:,l_fs_idx{kk},kk),m_fs_num(kk));
    m_fs_pos{kk} = des_prm.pos(m_fs_idx{kk},:);

    %% Draw figures

    if ismember(prm.freq(kk), figsave_freq) == 1
        fig_pos = [0, 0, 300, 350];

        %Reg
        fig(41)=figure(41);
        set(fig(41),'Position',fig_pos);
        hold on;
        plot(m_reg_pos{kk}(:,1),m_reg_pos{kk}(:,2),'xk','MarkerSize',4); %Microphones
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

        %Rand
        fig(42)=figure(42);
        set(fig(42),'Position',fig_pos);
        hold on;
        plot(m_rnd_pos{kk}(:,1),m_rnd_pos{kk}(:,2),'xk','MarkerSize',4); %Microphones
        plot(lc_prm.pos(:,1),lc_prm.pos(:,2),'.k','MarkerSize',1);
        plot(l_rnd_pos{kk}(:,1),l_rnd_pos{kk}(:,2),'ok','MarkerSize',4,'MarkerFaceColor','k'); %Loudspeakers
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

        %GSO
        fig(43)=figure(43);
        set(fig(43),'Position',fig_pos);
        hold on;
        plot(m_gso_pos{kk}(:,1),m_gso_pos{kk}(:,2),'xk','MarkerSize',4); %Microphones
        plot(lc_prm.pos(:,1),lc_prm.pos(:,2),'.k','MarkerSize',1);
        plot(l_gso_pos{kk}(:,1),l_gso_pos{kk}(:,2),'ok','MarkerSize',4,'MarkerFaceColor','k'); %Loudspeakers
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

        %Det
        fig(44)=figure(44);
        set(fig(44),'Position',fig_pos);
        hold on;
        plot(m_det_pos{kk}(:,1),m_det_pos{kk}(:,2),'xk','MarkerSize',4); %Microphones
        plot(lc_prm.pos(:,1),lc_prm.pos(:,2),'.k','MarkerSize',1);
        plot(l_det_pos{kk}(:,1),l_det_pos{kk}(:,2),'ok','MarkerSize',4,'MarkerFaceColor','k'); %Loudspeakers
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

        %MI
        fig(45)=figure(45);
        set(fig(45),'Position',fig_pos);
        hold on;
        plot(m_mi_pos{kk}(:,1),m_mi_pos{kk}(:,2),'xk','MarkerSize',4); %Microphones
        plot(lc_prm.pos(:,1),lc_prm.pos(:,2),'.k','MarkerSize',1);
        plot(l_mi_pos{kk}(:,1),l_mi_pos{kk}(:,2),'ok','MarkerSize',4,'MarkerFaceColor','k'); %Loudspeakers
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

        %FS
        fig(46)=figure(46);
        set(fig(46),'Position',fig_pos);
        hold on;
        plot(m_fs_pos{kk}(:,1),m_fs_pos{kk}(:,2),'xk','MarkerSize',4); %Microphones
        plot(lc_prm.pos(:,1),lc_prm.pos(:,2),'.k','MarkerSize',1);
        plot(l_fs_pos{kk}(:,1),l_fs_pos{kk}(:,2),'ok','MarkerSize',4,'MarkerFaceColor','k'); %Loudspeakers
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

        %EIM
        fig(47)=figure(47);
        set(fig(47),'Position',fig_pos);
        hold on;
        plot(m_eim_pos{kk}(:,1),m_eim_pos{kk}(:,2),'xk','MarkerSize',4); %Microphones
        plot(lc_prm.pos(:,1),lc_prm.pos(:,2),'.k','MarkerSize',1);
        plot(l_eim_pos{kk}(:,1),l_eim_pos{kk}(:,2),'ok','MarkerSize',4,'MarkerFaceColor','k'); %Loudspeakers
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
            fname = cell(7,1);

            fname{1} = sprintf('%s/fig6a',dir_results);
            fname{2} = sprintf('%s/fig6b',dir_results);
            fname{3} = sprintf('%s/fig6c',dir_results);
            fname{4} = sprintf('%s/fig6d',dir_results);
            fname{5} = sprintf('%s/fig6e',dir_results);
            fname{6} = sprintf('%s/fig6f',dir_results);
            fname{7} = sprintf('%s/fig6g',dir_results);
            save_figures(fig(41:47), fname, fig_type);
        end
    end
end
