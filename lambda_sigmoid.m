%pre-setting figure properties
fn='Helvetica';fs_labels=16;fs_axis=16;wd=20;ht=18;

%A_star range
A_star = linspace(0,1,1000);
% lambda(A_star) function
l = @(A_star) 1./(1+exp(-(16*A_star-8)));

%figure plot
sigmoid_fig=figure('Name','Sigmoidal Lambda');clf;
box on;hold on;grid on;
xlim([0 1]);ylim([0 1]);
plot(A_star,l(A_star),'m-','LineWidth',2);
%x1-axis settings
ax1=gca;
% % ax1.XLabel.String='$A^{*}$';
% % ax1.XLabel.Interpreter='latex';
% % ax1.XLabel.FontSize = fs_labels;
ax1.XTick = 0:0.1:1;
% % ax1.XLabel.FontSize = fs_axis;
% % ax1.XLabel.FontName = fn;
% % ax1.YLabel.String='$\lambda(A^{*})$';
% % ax1.YLabel.Interpreter='latex';
ax1.YTick = 0:0.1:1;
% % ax1.YLabel.FontSize = fs_axis;
% % ax1.YLabel.FontName = fn;
ax1.FontSize = fs_labels;
ax1.TickDir = 'out';

ax1_pos = ax1.Position;

ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
ax2.YTick=[];
ax2.XTick = 0:0.2:1;

valueOfA = zeros(1,6);

for i = 1:6
    val = ax2.XTick(i)*2760;
    valueOfA(i) = val;
end
ax2.XTickLabel = {valueOfA(1),valueOfA(2),valueOfA(3),valueOfA(4),valueOfA(5),valueOfA(6)};
% {valueOfA(1),valueOfA(2),valueOfA(3),valueOfA(4),valueOfA(5),valueOfA(6),valueOfA(7),valueOfA(8),valueOfA(9),valueOfA(10),valueOfA(11)};
ax2.TickDir = 'out';
ax2.FontSize = fs_axis;

set(gcf,'Units','centimeters','Position',[0 0 wd ht],'PaperUnits','centimeters','PaperSize',[wd+2 ht+2]);
hold off;

%saving figure
epsFileName = 'figures\\sigmoidal_lambda.eps';fullFileName =fullfile(epsFileName);print(sigmoid_fig,fullFileName,'-depsc');
tiffFileName = 'figures\\sigmoidal_lambda.tiff';fullFileName2 = fullfile(tiffFileName);print(sigmoid_fig,fullFileName2,'-dtiff');

