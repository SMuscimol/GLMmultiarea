function out = plotweights_exp(n, pMax, alpha, with_mean, mean_subtracted, varargin)
    
    if nargin < 5
        mean_subtracted = false;
    end
    
    if nargin < 4
        with_mean = false;
    end
    
    if nargin < 3
        alpha=1;
    end
    
    if nargin < 2
        pMax = 0.05;
    end

    hold on
    if with_mean
        m = exp(n.wPlain(1))/0.1;
    else
        m = 1;
    end
    
    if mean_subtracted
        subt = 1;
    else
        subt = 0.;
    end
    
    try
        str = '#4ba933';
        c= sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        c = cat(2, c, alpha);
        out.LED = m .* ([n.psWaldAv.LEDOnset.data; n.psWaldAv.LEDDelay.data]<pMax) .* ([exp(n.ws.LEDOnset.data); exp(n.ws.LEDDelay.data)]-subt);
        plot(0.05:0.1:0.95, out.LED, 'DisplayName','LED','color',c, varargin{:});
        
    catch
        warning('LED not found')
    end
     
    try
        str = '#ed801c';
        c= sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        c = cat(2, c, alpha);
        out.whisk =  m .* ([n.psWaldAv.wOnset.data; n.psWaldAv.wDelay.data]<pMax) .* ([exp(n.ws.wOnset.data); exp(n.ws.wDelay.data)]-subt);
        plot(1.05:0.1:1.95,out.whisk, 'DisplayName','whisker stim','color',c,varargin{:})
    catch
        warning('whisker onset not found')
     end
    
    try
        str = '#e22520';
        c= sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        c = cat(2, c, alpha);
        out.sound = m .* ([n.psWaldAv.soundOnset.data; n.psWaldAv.soundDelay.data]<pMax) .* ([exp(n.ws.soundOnset.data); exp(n.ws.soundDelay.data)]-subt);
        plot(2.05:0.1:2.95, out.sound, 'DisplayName','sound','color',c,varargin{:})
    catch
        warning('sound not found')
    end
   
    try
        c = [0.1,0.1,0.1];
        c = cat(2, c, alpha);
        out.lickOnset1 = m .* ([n.psWaldAv.lickOnsetPre.data; n.psWaldAv.lickOnsetPost.data(1)]<pMax) .* ([exp(n.ws.lickOnsetPre.data);exp(n.ws.lickOnsetPost.data(1))]-subt);
        out.lickOnset2 = m .* ([n.psWaldAv.lickOnsetPost.data; n.psWaldAv.lickOnsetDelay.data]<pMax) .* ([exp(n.ws.lickOnsetPost.data); exp(n.ws.lickOnsetDelay.data)]-subt);
        plot(2.05:0.1:2.25, out.lickOnset1, 'color',c,varargin{:},'LineStyle','--')
        plot(2.25:0.1:3.15, out.lickOnset2, 'DisplayName','lick onset','color',c,varargin{:})
    catch
        warning('lick onset not found')
    end
    xlabel('t [s]');
    ylabel('exp(X\beta)');
    
end

