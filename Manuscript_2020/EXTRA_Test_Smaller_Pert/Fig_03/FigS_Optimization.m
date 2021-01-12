%% Figure S5. Optimizing a controller for different subsystems.
%   Gómez-Schiavon & El-Samad
%   July 2020
clear

%% Panel B
myXL = 'Y synthesis rate (\mu_Y)';
myYL = 'CoRa_{\mu_Y\in\Theta}(\mu_Y)';

motifs = {'FADv1','FADv1p2','FADv1p3'};
for m = 1:length(motifs)
    load(cat(2,'DATA_DYms_',motifs{m},'_FigS5_mY_mY.mat'))
    pp = fieldnames(p);
    pN = regexprep(pp,'^m','\\mu_'); pN = regexprep(pN,'^g','\\gamma_'); pN = regexprep(pN,'^e','\\eta_'); pN = regexprep(pN,'_P','_+'); pN = regexprep(pN,'_M','_-');pN = regexprep(pN,'^b','\\beta_');pN = regexprep(pN,'kD','K_D');pN = regexprep(pN,'_$','');pN = regexprep(pN,'_Us','_\{Us\}');
    rhoN = regexprep(rho.name,'^m','\\mu_'); rhoN = regexprep(rhoN,'^g','\\gamma_'); rhoN = regexprep(rhoN,'^e','\\eta_'); rhoN = regexprep(rhoN,'_P','_+'); rhoN = regexprep(rhoN,'_M','_-');rhoN = regexprep(rhoN,'^b','\\beta_');rhoN = regexprep(rhoN,'kD','K_D');rhoN = regexprep(rhoN,'_$','');rhoN = regexprep(rhoN,'_Us','_\{Us\}');

    for i = 1:2
        fig = figure();
        fig.Units = 'inches';
        fig.PaperPosition = [2 1 4 2];
        fig.Position = fig.PaperPosition;
        C = [0 1 1;0 0.962745070457458 1;0 0.925490200519562 1;0 0.888235330581665 1;0 0.850980401039124 1;0 0.813725471496582 1;0 0.776470601558685 1;0 0.739215731620789 1;0 0.701960802078247 1;0 0.664705872535706 1;0 0.627451002597809 1;0 0.596078455448151 1;0 0.564705908298492 1;0 0.533333361148834 1;0 0.501960813999176 1;0 0.470588266849518 1;0 0.439215689897537 1;0 0.407843142747879 1;0 0.376470595598221 1;0 0.345098048448563 1;0 0.313725501298904 1;0 0.285205006599426 0.909090936183929;0 0.256684511899948 0.818181812763214;0 0.228164002299309 0.727272748947144;0 0.199643507599831 0.636363625526428;0 0.171122997999191 0.545454561710358;0 0.142602503299713 0.454545468091965;0 0.114082001149654 0.363636374473572;0 0.0855614989995956 0.272727280855179;0 0.0570410005748272 0.181818187236786;0 0.0285205002874136 0.0909090936183929;0 0 0;0 0 0;0.0909090936183929 0.0142602501437068 0.0909090936183929;0.181818187236786 0.0285205002874136 0.181818187236786;0.272727280855179 0.0427807494997978 0.272727280855179;0.363636374473572 0.0570410005748272 0.363636374473572;0.454545468091965 0.0713012516498566 0.454545468091965;0.545454561710358 0.0855614989995956 0.545454561710358;0.636363625526428 0.0998217537999153 0.636363625526428;0.727272748947144 0.114082001149654 0.727272748947144;0.818181812763214 0.128342255949974 0.818181812763214;0.909090936183929 0.142602503299713 0.909090936183929;1 0.156862750649452 1;1 0.191721141338348 1;1 0.226579532027245 1;1 0.261437922716141 1;1 0.296296298503876 1;1 0.331154674291611 1;1 0.366013079881668 1;1 0.400871455669403 1;1 0.43572986125946 1;1 0.470588237047195 1;1 0.499108731746674 1;1 0.527629256248474 1;1 0.55614972114563 1;1 0.58467024564743 1;1 0.613190710544586 1;1 0.641711235046387 1;1 0.670231759548187 1;1 0.698752224445343 1;1 0.727272748947144 1;1 0.755793213844299 1;1 0.7843137383461 1];
        hold on;
            [a b] = unique(p.(pp{i}));
            I = [b(1):b(end)];
            for j = 1:length(I)
                plot(rho.values,DYs(I(j),:),...
                    'LineWidth',3,'Color',C(1+((j-1)*10),:),...
                    'DisplayName',cat(2,pN{i},'_0\times10^{',num2str(log10(a(j)/a(ceil(length(I)/2)))),'}'))
            end
                ylabel(myYL)
                ylim([0 1])
                xlabel(myXL)
                xlim([min(rho.values) max(rho.values)])
                title(cat(2,'Changing ',pN{i}),'FontSize',12)
                set(gca,'XScale','log','XTick',10.^[-2:2:2],'FontSize',12)
                box on

                annotation('textbox',[0.6 0.3 0.25 0.125],'BackgroundColor',[1 1 1],...
                    'String',{cat(2,pN{i},'= ',num2str(mode(p.(pp{i}))))},...
                    'FitBoxToText','on','LineWidth',3,'FontSize',10);
            print(gcf,cat(2,'RAW_FigS5_',motifs{m},'_',pp{i}),'-dpng','-r300')
    end
end

%% Panel C
myL = 'CoRa_{\mu_Y\in\Theta}(\mu_Y)';
motifs = {'FADv1','FADv1p2','FADv1p3'};
for m = 1:length(motifs)
    load(cat(2,'DATA_OptDY_',motifs{m},'_FigS5_mY_mY.mat'))
    load('DATA_DYms_FADv1_Fig3_mY_mY.mat','rho')    % Just for rho values
    % Optimization
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 2];
    fig.Position = fig.PaperPosition;
    set(fig,'defaultAxesColorOrder',[0 0 0; 0.6 0.6 0.6]);
    hold on;
        plot(10.^DYs.thr,'LineWidth',3,...
            'LineWidth',3,'Color',[0 0 0])
            ylabel(cat(2,'|CoRa\leq',epsT,'|'),'FontSize',12)
            ylim(10.^[0 6])
            set(gca,'YScale','log',...
                'YTick',10.^[-6:2:6],'FontSize',12)
    yyaxis right
        hold on
        plot(DYs.min,'LineWidth',3,...
            'LineWidth',3,'Color',[0.6 0.6 0.6])
            ylabel(cat(2,'min(CoRa)'),'FontSize',12)
            xlabel('Iteration')
            ylim(10.^[-5 1])
            xlim([0 1000])
            box on
            set(gca,'YScale','log','XGrid','on','XTick',[0:250:1000],...
                'YTick',10.^[-6:2:6],'FontSize',12)
        print(gcf,cat(2,'RAW_FigS5_',motifs{m},'_Opt'),'-dpng','-r300')
    % Parameters
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 2];
    fig.Position = fig.PaperPosition;
    hold on;
        plot(ps.mU,'LineWidth',3,...
            'DisplayName','\mu_U')
        plot(ps.mW,'LineWidth',3,...
            'DisplayName','\mu_W')
        plot(ps.eP,'LineWidth',3,...
            'DisplayName','\eta_+')
            ylim(10.^[-3 3.1])
            xlabel('Iteration','FontSize',12)
            ylabel('Parameter','FontSize',12)
            xlim([0 1000])
            box on
            set(gca,'YScale','log','XGrid','on','XTick',[0:250:1000],...
                'YTick',10.^[-6:2:6],'FontSize',12)
            legend('Location','NorthWest')
        print(gcf,cat(2,'RAW_FigS5_',motifs{m},'_Par'),'-dpng','-r300')
    % CoRa lines
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 4 2];
    fig.Position = fig.PaperPosition;
    C = [0.619607865810394 0.0235294122248888 0.380392163991928;0.622471213340759 0.0358543433248997 0.387052595615387;0.625334620475769 0.0481792725622654 0.393713057041168;0.628197968006134 0.0605042017996311 0.400373488664627;0.631061315536499 0.0728291347622871 0.407033920288086;0.633924722671509 0.0851540639996529 0.413694381713867;0.636788070201874 0.0974789932370186 0.420354813337326;0.639651417732239 0.109803922474384 0.427015244960785;0.642514824867249 0.12212885171175 0.433675706386566;0.645378172397614 0.134453788399696 0.440336138010025;0.648241519927979 0.146778717637062 0.446996569633484;0.651104927062988 0.159103646874428 0.453657031059265;0.653968274593353 0.171428576111794 0.460317462682724;0.656831622123718 0.183753505349159 0.466977924108505;0.659695029258728 0.196078434586525 0.473638355731964;0.662558376789093 0.208403363823891 0.480298787355423;0.665421724319458 0.220728293061256 0.486959248781204;0.668285131454468 0.233053222298622 0.493619680404663;0.671148478984833 0.245378151535988 0.500280141830444;0.674011826515198 0.257703095674515 0.506940543651581;0.676875233650208 0.27002802491188 0.513601005077362;0.679738581180573 0.282352954149246 0.520261466503143;0.682601928710938 0.294677883386612 0.52692186832428;0.685465335845947 0.307002812623978 0.533582329750061;0.688328683376312 0.319327741861343 0.540242791175842;0.691192030906677 0.331652671098709 0.546903192996979;0.694055438041687 0.343977600336075 0.55356365442276;0.696918785572052 0.356302529573441 0.560224115848541;0.699782133102417 0.368627458810806 0.566884517669678;0.702645540237427 0.380952388048172 0.573544979095459;0.705508887767792 0.393277317285538 0.58020544052124;0.708372235298157 0.405602246522903 0.586865842342377;0.711235642433167 0.417927175760269 0.593526303768158;0.714098989963531 0.430252104997635 0.600186765193939;0.716962337493896 0.442577034235001 0.606847167015076;0.719825744628906 0.454901963472366 0.613507628440857;0.722689092159271 0.467226892709732 0.620168089866638;0.725552439689636 0.479551821947098 0.626828491687775;0.728415846824646 0.491876751184464 0.633488953113556;0.731279194355011 0.504201710224152 0.640149414539337;0.734142541885376 0.516526639461517 0.646809816360474;0.737005949020386 0.528851568698883 0.653470277786255;0.739869296550751 0.541176497936249 0.660130739212036;0.742732644081116 0.553501427173615 0.666791200637817;0.745596051216125 0.56582635641098 0.673451602458954;0.74845939874649 0.578151285648346 0.680112063884735;0.751322746276855 0.590476214885712 0.686772525310516;0.754186153411865 0.602801144123077 0.693432927131653;0.75704950094223 0.615126073360443 0.700093388557434;0.759912848472595 0.627451002597809 0.706753849983215;0.762776255607605 0.639775931835175 0.713414251804352;0.76563960313797 0.65210086107254 0.720074713230133;0.768502950668335 0.664425790309906 0.726735174655914;0.771366357803345 0.676750719547272 0.733395576477051;0.77422970533371 0.689075648784637 0.740056037902832;0.777093052864075 0.701400578022003 0.746716499328613;0.779956459999084 0.713725507259369 0.75337690114975;0.782819807529449 0.726050436496735 0.760037362575531;0.785683155059814 0.7383753657341 0.766697824001312;0.788546562194824 0.750700294971466 0.773358225822449;0.791409909725189 0.763025224208832 0.78001868724823;0.794273257255554 0.775350153446198 0.786679148674011;0.797136664390564 0.787675082683563 0.793339550495148;0.800000011920929 0.800000011920929 0.800000011920929]; C = C(end:-1:1,:);
    hold on;
    for i = [1:200:size(DYs.curv,1)]
        i = round(i);
        plot(rho.values,DYs.curv(i,:),'LineWidth',3,'Color',C(ceil(i/16),:),...
            'DisplayName',cat(2,'Iteration #',num2str(i-1)))
    end
        plot(rho.values,DYs.curv(end,:),'LineWidth',3,'Color',[0 0 0],...
            'DisplayName',cat(2,'Iteration #',num2str(size(DYs.curv,1)-1)))

            ylabel(myL,'FontSize',12)
            ylim([0 1])
            xlabel('Y synthesis rate (\mu_Y)','FontSize',12)
            xlim([0.001 1000])
            box on
            set(gca,'XScale','log','XTick',10.^[-2:2:2],'FontSize',12)
        print(gcf,cat(2,'RAW_FigS5_',motifs{m},'_CoRa'),'-dpng','-r300')
end
    