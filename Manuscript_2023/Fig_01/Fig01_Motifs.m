%% Figure 1. CoRa provides a unifying framework to compare different feedback control motifs
%   Gómez-Schiavon & El-Samad
%   March 2023
clear

%% All panels
motifs = {'ATFv1','ATFv2','FADv1','FADv2','BNFv1','BNFv2','FDPv1','FDPv2','FFLv1','BMFv1','BMFv2'};
for m = 1:length(motifs)
    load(cat(2,'DATA_CoRams_',motifs{m},'_Fig1_mY_mY.mat'))

    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 1 3 3];
    fig.Position = fig.PaperPosition;
    C = [0 1 1;0 0.962745070457458 1;0 0.925490200519562 1;0 0.888235330581665 1;0 0.850980401039124 1;0 0.813725471496582 1;0 0.776470601558685 1;0 0.739215731620789 1;0 0.701960802078247 1;0 0.664705872535706 1;0 0.627451002597809 1;0 0.596078455448151 1;0 0.564705908298492 1;0 0.533333361148834 1;0 0.501960813999176 1;0 0.470588266849518 1;0 0.439215689897537 1;0 0.407843142747879 1;0 0.376470595598221 1;0 0.345098048448563 1;0 0.313725501298904 1;0 0.285205006599426 0.909090936183929;0 0.256684511899948 0.818181812763214;0 0.228164002299309 0.727272748947144;0 0.199643507599831 0.636363625526428;0 0.171122997999191 0.545454561710358;0 0.142602503299713 0.454545468091965;0 0.114082001149654 0.363636374473572;0 0.0855614989995956 0.272727280855179;0 0.0570410005748272 0.181818187236786;0 0.0285205002874136 0.0909090936183929;0 0 0;0 0 0;0.0909090936183929 0.0142602501437068 0.0909090936183929;0.181818187236786 0.0285205002874136 0.181818187236786;0.272727280855179 0.0427807494997978 0.272727280855179;0.363636374473572 0.0570410005748272 0.363636374473572;0.454545468091965 0.0713012516498566 0.454545468091965;0.545454561710358 0.0855614989995956 0.545454561710358;0.636363625526428 0.0998217537999153 0.636363625526428;0.727272748947144 0.114082001149654 0.727272748947144;0.818181812763214 0.128342255949974 0.818181812763214;0.909090936183929 0.142602503299713 0.909090936183929;1 0.156862750649452 1;1 0.191721141338348 1;1 0.226579532027245 1;1 0.261437922716141 1;1 0.296296298503876 1;1 0.331154674291611 1;1 0.366013079881668 1;1 0.400871455669403 1;1 0.43572986125946 1;1 0.470588237047195 1;1 0.499108731746674 1;1 0.527629256248474 1;1 0.55614972114563 1;1 0.58467024564743 1;1 0.613190710544586 1;1 0.641711235046387 1;1 0.670231759548187 1;1 0.698752224445343 1;1 0.727272748947144 1;1 0.755793213844299 1;1 0.7843137383461 1];
    hold on;
        pp = fieldnames(p); pp = pp{1};
        pN = regexprep(pp,'^m','\\mu_'); pN = regexprep(pN,'^g','\\gamma_'); pN = regexprep(pN,'^e','\\eta_'); pN = regexprep(pN,'_P','_+'); pN = regexprep(pN,'_M','_-');pN = regexprep(pN,'^b','\\beta_');pN = regexprep(pN,'kD','K_D');pN = regexprep(pN,'_$','');pN = regexprep(pN,'_Us','_\{Us\}');
        rhoN = regexprep(rho.name,'^m','\\mu_'); rhoN = regexprep(rhoN,'^g','\\gamma_'); rhoN = regexprep(rhoN,'^e','\\eta_'); rhoN = regexprep(rhoN,'_P','_+'); rhoN = regexprep(rhoN,'_M','_-');rhoN = regexprep(rhoN,'^b','\\beta_');rhoN = regexprep(rhoN,'kD','K_D');rhoN = regexprep(rhoN,'_$','');rhoN = regexprep(rhoN,'_Us','_\{Us\}');
        for j = 1:7
            plot(rho.values,CoRas(j,:),'DisplayName','CoRa(\mu_Y)',...
                'LineWidth',3,'Color',C(1+((j-1)*10),:))
        end
            xlabel('Y synthesis rate (\mu_Y)','FontSize',12)
            xlim([0.001 1000])
            ylabel('CoRa_{\mu_Y\in\Theta}(\mu_Y)','FontSize',12)
            ylim([0 1])
            if(strcmp(pp,'mW'))
                title('Changing W synthesis rate','FontSize',12)
            elseif(strcmp(pp,'mU'))
                title('Changing U synthesis rate','FontSize',12)
            end
            set(gca,'XScale','log','XTick',10.^[-2:2:2],'FontSize',12)
            box on

            annotation('textbox',[0.6 0.3 0.25 0.125],...
                'String',{cat(2,pN,'= ',num2str(p.(pp)(4)))},...
                'FitBoxToText','on','LineWidth',3,'FontSize',12);
        print(gcf,cat(2,'RAW_Fig1_',motifs{m}),'-dpng','-r300')
end