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
pos_rndi_gen_flg = 0; % Flag for regenerating random sampling of internal points

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

%% Sampling inside internal region

m_int_num = 6; %Number of internal control points
m_eimi_num_max = 10; %Maximum number of internal control points for Reg+EIMi
m_rndi_ntrial = 1000000; %Number of trials of random sampling

%Indexes of interior region of target region
des_int_idx = find(des_prm.pos(:,1)<0.38 & des_prm.pos(:,1)>-0.38 & des_prm.pos(:,2)<0.48 & des_prm.pos(:,2)>-0.48);
des_int_num = length(des_int_idx);

%Files of random sampling
if rev_prm.flg == 0
    dirname_rndi_pos = sprintf('../data/pos_int/freefield_x%d_y%d/int%d_trial%d',des_prm.len_x*100, des_prm.len_y*100,m_int_num,m_rndi_ntrial);
elseif rev_prm.flg == 1
    dirname_rndi_pos = sprintf('../data/pos_int/abs%.2f_x%d_y%d/int%d_trial%d',rev_prm.abs_ratio,des_prm.len_x*100, des_prm.len_y*100,m_int_num,m_rndi_ntrial);
end

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

%Reg+CHIEF
m_regchief_num = zeros(prm.freq_len,1);
m_regchief_idx = cell(prm.freq_len,1);
m_regchief_pos = cell(prm.freq_len,1);

%Reg+EIM
m_regeim_num = zeros(prm.freq_len,1);
m_regeim_idx = cell(prm.freq_len,1);
m_regeim_pos = cell(prm.freq_len,1);

%Reg+EIMi
m_regeimin_num = zeros(prm.freq_len,m_eimi_num_max);
m_regeimin_idx = cell(prm.freq_len,m_eimi_num_max);
m_regeimin_pos = cell(prm.freq_len,m_eimi_num_max);

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

    %% Regular - EIM (Reg+EIM)
    
    %Microphones
    m_regeim_num(kk) = m_eim_num(kk);
    [~, m_regeim_idx{kk}] = lmp_s_eim_n(H_lc2d(:, l_reg_idx{kk},kk), m_regeim_num(kk));
    m_regeim_pos{kk} = des_prm.pos(m_regeim_idx{kk},:);

    %% Regular - Regular + EIM for additional points (Reg+EIMi)

    %Microphones
    for jj=1:m_eimi_num_max
        m_regeimin_num(kk,jj) = m_eim_num(kk);
        m_eimi_num = jj;
        if m_eimi_num > m_regeimin_num(kk,jj)
            [~, m_eimin_idx] = lmp_s_eim_n(H_lc2d(:, l_reg_idx{kk},kk), m_regeimin_num(kk));
            m_regeimin_idx{kk,jj} = m_eimin_idx;
            m_regeimin_pos{kk,jj} = des_prm.pos(m_regeimin_idx{kk,jj},:);
        else
            m_regin_num = m_regeimin_num(kk,jj)-m_eimi_num;
            m_regin_pos_tmp = rect_perim(m_regin_num, des_prm.len_x, des_prm.len_y, 0, 0);
            m_regin_idx = zeros(m_regin_num,1);
            for ii=1:m_regin_num
                dd = sqrt((m_regin_pos_tmp(ii,1)-des_prm.pos(:,1)).^2 + (m_regin_pos_tmp(ii,2)-des_prm.pos(:,2)).^2);
                [~, m_regin_idx(ii)] = min(dd); 
            end

            ker = null(H_lc2d(m_regin_idx, l_reg_idx{kk},kk));
            ker_w = H_lc2d(:, l_reg_idx{kk}, kk) * ker;
            [~, m_eimi_idx] = lmp_s_eim_n(ker_w, m_eimi_num);
            m_regeimin_idx{kk,jj} = [m_regin_idx; m_eimi_idx];
            m_regeimin_pos{kk,jj} = des_prm.pos(m_regeimin_idx{kk,jj},:);
        end
    end

    %% Regular - Regular + Random sampling for additional points

    fname_rndi_pos = sprintf('%s/pos_rndi_f%d.mat',dirname_rndi_pos,prm.freq(kk));

    if pos_rndi_gen_flg == 1
        m_regrndi_num = m_eim_num(kk);
        m_rndi_num = m_int_num;
        m_regrndi_idx = zeros(m_regrndi_num,m_rndi_ntrial);
        m_regrndi_pos = zeros(m_regrndi_num,2,m_rndi_ntrial);

        if m_rndi_num>m_regrndi_num
            for ii=1:m_rndi_ntrial
                m_rndi_idx = des_int_idx(randsample(des_prm.num,m_rndi_num));
                m_regrndi_idx(:,ii) = m_rndi_idx;
                m_regrndi_pos(:,:,ii) = des_prm.pos(m_regrndi_idx(:,ii),:);
            end
        else
            m_regin_num = m_regrndi_num-m_int_num;
            m_regin_pos_tmp = rect_perim(m_regin_num, des_prm.len_x, des_prm.len_y, 0, 0);
            m_regin_idx = zeros(m_regin_num,1);
            for ii=1:m_regin_num
                dd = sqrt((m_regin_pos_tmp(ii,1)-des_prm.pos(:,1)).^2 + (m_regin_pos_tmp(ii,2)-des_prm.pos(:,2)).^2);
                [~, m_regin_idx(ii)] = min(dd); 
            end

            des_regi_idx = (1:des_prm.num)';
            des_regi_idx(m_regin_idx) = [];
            des_regi_num = length(des_regi_idx);

            for ii=1:m_rndi_ntrial
                m_rndi_idx = des_regi_idx(randsample(des_regi_num,m_rndi_num));
                m_regrndi_idx(:,ii) = [m_regin_idx; m_rndi_idx];
                m_regrndi_pos(:,:,ii) = des_prm.pos(m_regrndi_idx(:,ii),:);
            end
        end

        save(fname_rndi_pos,'m_regrndi_num','m_regrndi_idx','m_regrndi_pos');
    elseif pos_rndi_gen_flg == 0
        load(fname_rndi_pos,'m_regrndi_num','m_regrndi_idx','m_regrndi_pos');
    end

    %% Regular - Regular + CHIEF (Reg+CHIEF)

    %Microphones
    m_regchief_num(kk) =  m_eim_num(kk);
    m_chief_num = 6;

    m_regin_num = m_regchief_num(kk)-m_chief_num;
    m_regin_pos_tmp = rect_perim(m_regin_num, des_prm.len_x, des_prm.len_y, 0, 0);
    m_regin_idx = zeros(m_regin_num,1);
    for ii=1:m_regin_num
        dd = sqrt((m_regin_pos_tmp(ii,1)-des_prm.pos(:,1)).^2 + (m_regin_pos_tmp(ii,2)-des_prm.pos(:,2)).^2);
        [~, m_regin_idx(ii)] = min(dd); 
    end

    m_chief6_c_pos_tmp = [(max(des_prm.pos(:,1))+min(des_prm.pos(:,1)))/2,(max(des_prm.pos(:,2))+min(des_prm.pos(:,2)))/2];
    dd = sqrt((m_chief6_c_pos_tmp(1)-des_prm.pos(:,1)).^2 + (m_chief6_c_pos_tmp(2)-des_prm.pos(:,2)).^2);
    [~, m_chief6_c_idx] = min(dd); 
    m_chief6_c_pos = des_prm.pos(m_chief6_c_idx,:);

    m_chief6_pos_tmp = zeros(6,2);
    m_chief6_pos_tmp(1,:) = m_chief6_c_pos;
    m_chief6_pos_tmp(2,:) = m_chief6_c_pos + [des_prm.dx,0];
    m_chief6_pos_tmp(3,:) = m_chief6_c_pos + [0,des_prm.dy*3];
    m_chief6_pos_tmp(4,:) = m_chief6_c_pos + [des_prm.dx,des_prm.dy*3];
    m_chief6_pos_tmp(5,:) = m_chief6_c_pos + [des_prm.dx*3,0];
    m_chief6_pos_tmp(6,:) = m_chief6_c_pos + [des_prm.dx*4,0];

    m_chief6_idx = zeros(6,1);
    for ii=1:6
        dd = sqrt((m_chief6_pos_tmp(ii,1)-des_prm.pos(:,1)).^2 + (m_chief6_pos_tmp(ii,2)-des_prm.pos(:,2)).^2);
        [~, m_chief6_idx(ii)] = min(dd);
    end
    m_regchief_idx{kk} = [m_regin_idx; m_chief6_idx];
    m_regchief_pos{kk} = des_prm.pos(m_regchief_idx{kk},:);

    %% Draw figures

    if ismember(prm.freq(kk), figsave_freq) == 1
        close all;
        fig_pos = [0, 0, 300, 350];

        %Reg
        fig(41)=figure(41);
        set(fig(41),'Position',fig_pos);
        hold on;
        plot(m_reg_pos{kk,1}(:,1),m_reg_pos{kk,1}(:,2),'xk','MarkerSize',4); %Microphones
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

        %Reg+CHIEF
        fig(43)=figure(43);
        set(fig(43),'Position',fig_pos);
        hold on;
        plot(m_regchief_pos{kk,1}(:,1),m_regchief_pos{kk,1}(:,2),'xk','MarkerSize',4); %Microphones
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

        %Reg+EIMi
        fig(45)=figure(45);
        set(fig(45),'Position',fig_pos);
        hold on;
        plot(m_regeimin_pos{kk,m_int_num}(:,1),m_regeimin_pos{kk,m_int_num}(:,2),'xk','MarkerSize',4); %Microphones
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

        %Reg+EIM
        fig(46)=figure(46);
        set(fig(46),'Position',fig_pos);
        hold on;
        plot(m_regeim_pos{kk}(:,1),m_regeim_pos{kk}(:,2),'xk','MarkerSize',4); %Microphones
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
            fname = cell(3,1);

            fname{1} = sprintf('%s/fig11a',dir_results);
            fname{2} = sprintf('%s/fig11b',dir_results);
            fname{3} = sprintf('%s/fig14',dir_results);
            save_figures([fig(43),fig(45:46)], fname, fig_type);
        end
    end

    clear 'm_regrndi_num' 'm_regrndi_idx' 'm_regrndi_pos';
end
