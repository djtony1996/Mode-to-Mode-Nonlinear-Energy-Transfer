clear,clc

load('M_180_avg_final.mat','M_4d_avg_final');

load('full180_mean.mat','dkx','dky','nx','ny');

axis_range = [2 4.5 1.2 3.5];

max_kx_index_path = 21;
max_ky_index_path = 21;

positive_max_kx_index = zeros(max_ky_index_path,max_kx_index_path);
positive_max_ky_index = zeros(max_ky_index_path,max_kx_index_path);
positive_max_magnitude = zeros(max_ky_index_path,max_kx_index_path);
negative_max_kx_index = zeros(max_ky_index_path,max_kx_index_path);
negative_max_ky_index = zeros(max_ky_index_path,max_kx_index_path);
negative_max_magnitude = zeros(max_ky_index_path,max_kx_index_path);

for kx_index = 1: max_kx_index_path
    for ky_index = 1: max_ky_index_path
        M_single = M_4d_avg_final(:,:,ky_index,kx_index);
        
        [temp_max, temp_index] = max(M_single .* (M_single>0),[],'all');
        [temp_ky_index,temp_kx_index] = ind2sub(size(M_single), temp_index);
        positive_max_kx_index(ky_index,kx_index) = temp_kx_index;
        positive_max_ky_index(ky_index,kx_index) = temp_ky_index;
        positive_max_magnitude(ky_index,kx_index) = temp_max; 

        [temp_min, temp_index] = min(M_single .* (M_single<0),[],'all');
        [temp_ky_index,temp_kx_index] = ind2sub(size(M_single), temp_index);
        negative_max_kx_index(ky_index,kx_index) = temp_kx_index;
        negative_max_ky_index(ky_index,kx_index) = temp_ky_index;
        negative_max_magnitude(ky_index,kx_index) = temp_min;

    end
end


kx_0 = 0.2;
ky_0 = 0.4;
% for plotting pre-multiplied spectra
kx_start = [0:dkx:(max_kx_index_path-1)*dkx]; kx_start(1) = kx_0;% for infinity
ky_start = [0:dky:(max_ky_index_path-1)*dky]; ky_start(1) = ky_0;% for infinity
lambda_x_start = 2*pi./kx_start;
lambda_y_start = 2*pi./ky_start;

[coor_x_start,coor_y_start] = meshgrid(lambda_x_start.*180,lambda_y_start.*180);
coor_x_start_log = log10(coor_x_start);
coor_y_start_log = log10(coor_y_start);

kx_end_positive = (positive_max_kx_index-1).*dkx;
kx_end_positive(kx_end_positive==0) = kx_0;
ky_end_positive = (positive_max_ky_index-1).*dky;
ky_end_positive(ky_end_positive==0) = ky_0;

coor_x_end_log_positive = 180*2*pi./kx_end_positive;
coor_y_end_log_positive = 180*2*pi./ky_end_positive;

kx_end_negative = (negative_max_kx_index-1).*dkx;
kx_end_negative(kx_end_negative==0) = kx_0;
ky_end_negative = (negative_max_ky_index-1).*dky;
ky_end_negative(ky_end_negative==0) = ky_0;

coor_x_end_log_negative = 180*2*pi./kx_end_negative;
coor_y_end_log_negative = 180*2*pi./ky_end_negative;

sin_positive = (coor_y_end_log_positive - coor_y_start_log)./sqrt((coor_x_end_log_positive - coor_x_start_log).^2 + (coor_y_end_log_positive - coor_y_start_log).^2);
cos_positive = (coor_x_end_log_positive - coor_x_start_log)./sqrt((coor_x_end_log_positive - coor_x_start_log).^2 + (coor_y_end_log_positive - coor_y_start_log).^2);

sin_negative = (coor_y_end_log_negative - coor_y_start_log)./sqrt((coor_x_end_log_negative - coor_x_start_log).^2 + (coor_y_end_log_negative - coor_y_start_log).^2);
cos_negative = (coor_x_end_log_negative - coor_x_start_log)./sqrt((coor_x_end_log_negative - coor_x_start_log).^2 + (coor_y_end_log_negative - coor_y_start_log).^2);

u_quiver_positive = positive_max_magnitude .* cos_positive;
v_quiver_positive = positive_max_magnitude .* sin_positive;

u_quiver_negative = negative_max_magnitude .* cos_negative;
v_quiver_negative = negative_max_magnitude .* sin_negative;

%%

FontSize = 20;
TickLength      = 0.03;
left_coordinate   = 0.10;
bottom_coordinate = 0.15;
plot_width        = 0.35;
plot_height       = 0.80;

tickarray_major   = log10([0.1 1 10 100 1000]);
tickarray_minor   = log10([0.2:0.1:0.9 2:9 20:10:90 200:100:900 2000:1000:9000]);

figure
set(gcf,'Position',[680         556        950         421])
ax1 = axes("Position",[left_coordinate bottom_coordinate plot_width plot_height]);
quiver(coor_x_start_log, coor_y_start_log, u_quiver_negative, v_quiver_negative,4,'Color','k','LineWidth',1)
set(gca,'xtick',tickarray_major)
set(gca,'xticklabel',{'$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'})
set(gca,'ytick',tickarray_major)
set(gca,'yticklabel',{'$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'})
xlabel('$\lambda_x^+$','Interpreter','Latex')
ylabel('$\lambda_y^+$','Interpreter','Latex')
set(gca,'YDir','normal')
set(gca,'Fontsize',FontSize)
ax_obj = gca;
ax_obj.TickLabelInterpreter = 'Latex';
ax_obj.TickLength = [TickLength,TickLength];
ax_obj.XMinorTick = 'on';
ax_obj.XAxis.MinorTickValues = tickarray_minor;
ax_obj.YMinorTick = 'on';
ax_obj.YAxis.MinorTickValues = tickarray_minor;
axis square
text(2.1,3.3,'(a)','FontSize',20,'Interpreter','Latex')
axis(axis_range)
text(4.4,1.1,'$\infty$','FontSize',20,'Interpreter','Latex')
text(1.75,3.5,'$\infty$','FontSize',20,'Interpreter','Latex')

% x2 = log10(fliplr(lambda_x_start));
% y2 = log10(fliplr(lambda_y_start));
% ax2 = axes("Position",[left_coordinate bottom_coordinate plot_width plot_height]);
% axis([x2(1) x2(end) y2(1) y2(end)])
% axis square
% xlabel('$\lambda_x$','Interpreter','Latex')
% ylabel('$\lambda_y$','Interpreter','Latex')
% set(gca,'xtick',tickarray_major)
% set(gca,'xticklabel',{'$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'})
% set(gca,'ytick',tickarray_major)
% set(gca,'yticklabel',{'$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'})
% set(gca,'Fontsize',20)
% ax_obj = gca;
% ax_obj.TickLabelInterpreter = 'Latex';
% ax_obj.TickLength = [TickLength,TickLength];
% ax_obj.XMinorTick = 'on';
% ax_obj.XAxis.MinorTickValues = tickarray_minor;
% ax_obj.YMinorTick = 'on';
% ax_obj.YAxis.MinorTickValues = tickarray_minor;
% ax2.XAxisLocation = 'top';
% ax2.YAxisLocation = 'right';
% ax2.Color = 'none';
% ax1.Box = 'off';
% ax2.Box = 'off';


%%
clear,clc

load('M_590_avg_final.mat','M_4d_avg_final');

load('full590_mean.mat','dkx','dky','nx','ny');

axis_range = [2 4.5 1.2 3.5];

max_kx_index_path = 51;
max_ky_index_path = 51;

positive_max_kx_index = zeros(max_ky_index_path,max_kx_index_path);
positive_max_ky_index = zeros(max_ky_index_path,max_kx_index_path);
positive_max_magnitude = zeros(max_ky_index_path,max_kx_index_path);
negative_max_kx_index = zeros(max_ky_index_path,max_kx_index_path);
negative_max_ky_index = zeros(max_ky_index_path,max_kx_index_path);
negative_max_magnitude = zeros(max_ky_index_path,max_kx_index_path);

for kx_index = 1: max_kx_index_path
    for ky_index = 1: max_ky_index_path
        M_single = M_4d_avg_final(:,:,ky_index,kx_index);
        
        [temp_max, temp_index] = max(M_single .* (M_single>0),[],'all');
        [temp_ky_index,temp_kx_index] = ind2sub(size(M_single), temp_index);
        positive_max_kx_index(ky_index,kx_index) = temp_kx_index;
        positive_max_ky_index(ky_index,kx_index) = temp_ky_index;
        positive_max_magnitude(ky_index,kx_index) = temp_max; 

        [temp_min, temp_index] = min(M_single .* (M_single<0),[],'all');
        [temp_ky_index,temp_kx_index] = ind2sub(size(M_single), temp_index);
        negative_max_kx_index(ky_index,kx_index) = temp_kx_index;
        negative_max_ky_index(ky_index,kx_index) = temp_ky_index;
        negative_max_magnitude(ky_index,kx_index) = temp_min;

    end
end


kx_0 = 0.2;
ky_0 = 0.4;
% for plotting pre-multiplied spectra
kx_start = [0:dkx:(max_kx_index_path-1)*dkx]; kx_start(1) = kx_0;% for infinity
ky_start = [0:dky:(max_ky_index_path-1)*dky]; ky_start(1) = ky_0;% for infinity
lambda_x_start = 2*pi./kx_start;
lambda_y_start = 2*pi./ky_start;

[coor_x_start,coor_y_start] = meshgrid(lambda_x_start.*590,lambda_y_start.*590);
coor_x_start_log = log10(coor_x_start);
coor_y_start_log = log10(coor_y_start);

kx_end_positive = (positive_max_kx_index-1).*dkx;
kx_end_positive(kx_end_positive==0) = kx_0;
ky_end_positive = (positive_max_ky_index-1).*dky;
ky_end_positive(ky_end_positive==0) = ky_0;

coor_x_end_log_positive = 590*2*pi./kx_end_positive;
coor_y_end_log_positive = 590*2*pi./ky_end_positive;

kx_end_negative = (negative_max_kx_index-1).*dkx;
kx_end_negative(kx_end_negative==0) = kx_0;
ky_end_negative = (negative_max_ky_index-1).*dky;
ky_end_negative(ky_end_negative==0) = ky_0;

coor_x_end_log_negative = 590*2*pi./kx_end_negative;
coor_y_end_log_negative = 590*2*pi./ky_end_negative;

sin_positive = (coor_y_end_log_positive - coor_y_start_log)./sqrt((coor_x_end_log_positive - coor_x_start_log).^2 + (coor_y_end_log_positive - coor_y_start_log).^2);
cos_positive = (coor_x_end_log_positive - coor_x_start_log)./sqrt((coor_x_end_log_positive - coor_x_start_log).^2 + (coor_y_end_log_positive - coor_y_start_log).^2);

sin_negative = (coor_y_end_log_negative - coor_y_start_log)./sqrt((coor_x_end_log_negative - coor_x_start_log).^2 + (coor_y_end_log_negative - coor_y_start_log).^2);
cos_negative = (coor_x_end_log_negative - coor_x_start_log)./sqrt((coor_x_end_log_negative - coor_x_start_log).^2 + (coor_y_end_log_negative - coor_y_start_log).^2);

u_quiver_positive = positive_max_magnitude .* cos_positive;
v_quiver_positive = positive_max_magnitude .* sin_positive;

u_quiver_negative = negative_max_magnitude .* cos_negative;
v_quiver_negative = negative_max_magnitude .* sin_negative;

%%
FontSize = 20;
TickLength      = 0.03;
left_coordinate1  = 0.55;
bottom_coordinate = 0.15;
plot_width        = 0.35;
plot_height       = 0.80;


tickarray_major   = log10([0.1 1 10 100 1000]);
tickarray_minor   = log10([0.2:0.1:0.9 2:9 20:10:90 200:100:900 2000:1000:9000]);

figure(1)
set(gcf,'Position',[680         556        950         421])
ax1 = axes("Position",[left_coordinate1 bottom_coordinate plot_width plot_height]);
quiver(coor_x_start_log, coor_y_start_log, u_quiver_negative, v_quiver_negative,8,'Color','k','LineWidth',1)
set(gca,'xtick',tickarray_major)
set(gca,'xticklabel',{'$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'})
set(gca,'ytick',tickarray_major)
set(gca,'yticklabel',{'$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'})
xlabel('$\lambda_x^+$','Interpreter','Latex')
ylabel('$\lambda_y^+$','Interpreter','Latex')
set(gca,'YDir','normal')
set(gca,'Fontsize',FontSize)
ax_obj = gca;
ax_obj.TickLabelInterpreter = 'Latex';
ax_obj.TickLength = [TickLength,TickLength];
ax_obj.XMinorTick = 'on';
ax_obj.XAxis.MinorTickValues = tickarray_minor;
ax_obj.YMinorTick = 'on';
ax_obj.YAxis.MinorTickValues = tickarray_minor;
axis square
text(2.1,3.3,'(b)','FontSize',20,'Interpreter','Latex')
axis(axis_range)
text(4.4,1.1,'$\infty$','FontSize',20,'Interpreter','Latex')
text(1.75,3.5,'$\infty$','FontSize',20,'Interpreter','Latex')

% x2 = log10(fliplr(lambda_x_start));
% y2 = log10(fliplr(lambda_y_start));
% ax2 = axes("Position",[left_coordinate1 bottom_coordinate plot_width plot_height]);
% axis([x2(1) x2(end) y2(1) y2(end)])
% axis square
% xlabel('$\lambda_x$','Interpreter','Latex')
% ylabel('$\lambda_y$','Interpreter','Latex')
% set(gca,'xtick',tickarray_major)
% set(gca,'xticklabel',{'$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'})
% set(gca,'ytick',tickarray_major)
% set(gca,'yticklabel',{'$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'})
% set(gca,'Fontsize',20)
% ax_obj = gca;
% ax_obj.TickLabelInterpreter = 'Latex';
% ax_obj.TickLength = [TickLength,TickLength];
% ax_obj.XMinorTick = 'on';
% ax_obj.XAxis.MinorTickValues = tickarray_minor;
% ax_obj.YMinorTick = 'on';
% ax_obj.YAxis.MinorTickValues = tickarray_minor;
% ax2.XAxisLocation = 'top';
% ax2.YAxisLocation = 'right';
% ax2.Color = 'none';
% ax1.Box = 'off';
% ax2.Box = 'off';
