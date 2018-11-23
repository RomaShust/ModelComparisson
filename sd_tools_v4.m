% set of functions for Sniff Dependent(sd)odor responses  project
% this analysis is published in "Sniff invariant odor coding"
% Shusterman et al. bioRxiv 174417; doi: https://doi.org/10.1101/174417


function sd = sd_tools_v4()
    
    sd.all_units       = @all_units;
    sd.binned_raster   = @binned_raster;
    sd.cell_classify   = @cell_classify;
    sd.create_unit     = @create_unit;
    sd.default_params  = @default_params;
    sd.file_names      = @file_names;
    
    sd.get_md          = @get_md;
    sd.get_doublets    = @get_doublets;
    sd.get_phase       = @get_phase;
    sd.get_sniff_index = @get_sniff_index; 
    sd.get_triplets    = @get_triplets;
    
    sd.has_response    = @has_response;
    sd.md_xcorr        = @md_xcorr;
    sd.md_xcorr_boot   = @md_xcorr_boot;
    sd.mean_warp_par   = @mean_warp_par;
    
    sd.pop_classify    = @pop_classify;
    sd.pop_vector      = @pop_vector;
    sd.pop_xcorr_boot  = @pop_xcorr_boot;
    
    sd.psth            = @psth;
    sd.raster          = @raster;
    sd.raster4regr     = @raster4regr;
    sd.regr            = @regr;
    sd.sharp_responses = @sharp_responses;
    sd.spCountFastVsSlow = @spCountFastVsSlow;
    sd.split_psth        = @split_psth;
    
    sd.version           = @version;
    sd.writeVersionInfo  = @writeVersionInfo;
    sd.warping           = @warping;
    sd.spikeAlignmentModel  = @spikeAlignmentModel;
    sd.bestL_Estimation     = @bestL_Estimation;
    sd.spikeRealignment     = @spikeRealignment;
    sd.timeShiftByVolume    = @timeShiftByVolume;
    sd.meanAndStd           = @meanAndStd;
    sd.plotRasterAndPSTH    = @plotRasterAndPSTH;
    
    sd.diagBoxPlot          = @diagBoxPlot;
    sd.diagHist             = @diagHist;
end

function ver = version()
% returns the version number and the date of last update
%                                           Roma            10 Feb 2016
    ver = 4;   
    disp('version: 4     date:   10 Dec 2016')
end

function params = default_params()
% returns default parameters
%                                           Roma            10 Feb 2016

    params.range       = [1,500];    % (msec) post inhalation range
    params.sniff_range = [100, 500];
    params.wnd         = 10;         % (msec) psth boxcar window
    params.step        = 4;          % (msec) psth resolution
    params.r_type      = 'rt';       % rw = warped, rt = real time
    params.part        = 0.3;        % part of the data that you want
    params.n_split     = 1;          % number of parts to split the data
    params.nt          = 920;        % length of data in ms from sniff onset 
                                     % for real time and inhalation warped cycles  
    params.correction  = 32;         % correction for filtered sniff signal

end

function fn = file_names(mouse, sess)
% returns the strucutre fn of file names for a given mouse and session
% janelia's data id in mouse_jdb folder
%                                           Roma            1 Mar 2016
global fn_disk

    if isempty(fn_disk)
        fn_disk = setup_disk();
    end

    if ismac
        fn.disk      = fullfile('/Volumes', fn_disk);
    else
        fn.disk      = fn_disk;
    end

    fn.db_main_fold    = fullfile(fn.disk, 'Experiments', 'mouse_jdb');
    fn.db_sess_fold    = fullfile(fn.db_main_fold, sprintf('m%03d_%02d', mouse, sess));
    fn.db_spikes       = fullfile(fn.db_sess_fold, sprintf('spikes_%03d_%02d.mat', mouse, sess));
    fn.db_rsm_data     = fullfile(fn.db_sess_fold, sprintf('rsm_data_%03d_%02d.mat', mouse, sess));
    fn.db_sess         = fullfile(fn.db_sess_fold, sprintf('sess_%03d_%02d.mat', mouse, sess));
    fn.sniff_waveforms = fullfile(fn.db_sess_fold, sprintf('sniff_%03d_%02d.mat', mouse, sess)); 

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function disk = setup_disk()
        if ismac
            list = dir('/Volumes');
            kd = 0; d = cell(0);
            for il = 1:length(list)
                if list(il).name(1) ~='.'
                    kd = kd + 1;
                    d{kd}= list(il).name;
                end
            end
        else
            import java.io.*; 
            f=File('');
            r=f.listRoots;
            d = cell(0);
            for i=1:numel(r)
                d{i} = char(r(i));
            end
        end
        disp(d)
 
        [sel, ~] = listdlg('Liststring', d, ...
            'SelectionMode',    'single', ...
            'ListSize',         [150, 50], ...
            'PromptString', 'Choose disk:');
        disk = d{sel}; 
        fprintf('disk: %s \n', disk)
    end

end

function md = get_md(unit)
% creates a structure of metadata from structure unit
%                                           Roma             12.28.2011

    md = struct('unit_ref',     0, ...          % unit number from structure unit        
                'mouse',        0, ...          % mouse number
                'sess',         0, ...          % session number
                'cell_num',     0, ...          % cell number id from spikes_mmm_ss.mat
                'valve',        0, ...          % valve number for a given odor presenation
                'conc',         0, ...          % concentration
                'stimID',       '', ...         % string: '%s_%06.2f', odor tag and concentration
                'ID',           '');            % meta data ID 'mmm_ss_ccc'

    
    id = 0;
    for iu = 1:length(unit)    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        u  = unit(iu);
        
        % set of different valves
        valve = setdiff(unique(u.valve),0);
        nv    = length(valve);

        for iv = 1:nv  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
            conc = setdiff(unique(u.conc(u.valve==valve(iv))),0);
            conc = sort(conc);
            nc   = length(conc);
            
%             if nc > 1,  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                for ic = 1:nc  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    id = id + 1;
                    fprintf('%d    unit: %d   valve: %d  conc: %06.2d \n', id, iu, valve(iv), conc(ic))
                    
                    md(id).unit_ref = iu;
                    md(id).mouse    = u.mouse;
                    md(id).sess     = u.sess;
                    md(id).cell_num = u.cell_num;
                    md(id).valve    = valve(iv);
                    md(id).conc     = conc(ic);
                    md(id).stimID   = sprintf('%s_%06.2f', u.odorID{valve(iv)}, conc(ic));
                    md(id).ID       = sprintf('%03d_%02d_%03d', u.mouse, u.sess, u.cell_num);
                end    % ic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%             end
        end    % iv ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end

end

function [ind, sn_dur, inh_dur] = get_sniff_index(unit, md, params, sniff)
% the sniff index and corrsponding sniff durations for a given md structure
%
% function returns the index of sniffs for a given md structure based on fields unit_ref,
% valve, conc and given sniff number, if sniff = 0 ind returns all background sniffs 
%
% sniff numbers are as follows:
% 0     : sniffs when no odor is presented
% 1-N   : sniffs when odor is on
% -1:-3 : sniffs immediately after odor is turned off
%
% specifying sniff 0.1, 0.2, etc, retrieves sniffs prior to odor onset


    if numel(md)>1; error('one md at a time pls'); end

    ku = md.unit_ref;
    
    if sniff>0 && sniff<1
        ind = (unit(ku).sniff_dur>=params.sniff_range(1))&...
            (unit(ku).sniff_dur<=params.sniff_range(2));
            
        pre = round(sniff*10);
        tmp = ind & (unit(ku).odorSniffNum == 1);
        tmp = find(tmp)-pre;
        ind(:) = false;
        ind(tmp) = true;
        
    elseif sniff>0
        ind =   (unit(ku).sniff_dur>=params.sniff_range(1))&...
                (unit(ku).sniff_dur<=params.sniff_range(2))&...
                (unit(ku).odorSniffNum == sniff)&...
                (unit(ku).valve == md.valve)&...
                (unit(ku).conc == double(md.conc));
    else
        ind =   (unit(ku).sniff_dur>=params.sniff_range(1))&...
                (unit(ku).sniff_dur<=params.sniff_range(2))&...
                (unit(ku).odorSniffNum == sniff);        
    end
    sn_dur = unit(ku).sniff_dur(ind);
    inh_dur = unit(ku).inh_dur(ind);
end

function [all,t_axis,h] = raster(unit, md, params, sniff, varargin)
% builds raster matrix and if handle to axes provided, 
% it will build raster plot
%
% [all,t_axis, g] = dc.raster(unit, md(1), params, 1, g, 'by_freq');
%
%                                           Roma       17 Feb 2016

    if length(md)>1; error('one md at a time pls'); end

    range       = params.range;    
    t_axis      = range(1):range(2);
    
    [ind, sn_dur] = get_sniff_index(unit, md, params, sniff);
    all = unit(md.unit_ref).(params.r_type)(t_axis,ind);
    
    if (nargin > 4)
        disp(varargin{1})
        if ishandle(varargin{1})
            if nargin == 6 && strcmp(varargin{2}, 'by_freq')
                [~, sn_ind] = sort(sn_dur);
                [tr,sp] = find(all(:, sn_ind));
            else
                [tr,sp] = find(all);
            end
            subplot(varargin{1});
            h = plot(tr+range(1)-1,sp,'.');
        end
    end
    
end

function [unit, sp_pos, tr_num, sn_dur, sn_ind, inh_dur, inh_ind] = raster4regr(unit, md, params, sniff, snWF, varargin)
% builds raster matrix and if handle to axes provided, 
% it will build raster plot
%
% [all,t_axis, g] = dc.raster(unit, md(1), params, 1);
%
%                                           Roma       26 Feb 2016
    if length(md)>1;  error('one md at a time pls');  end

    range       = params.range;    
    t_axis      = range(1):range(2);
    
    [ind, sn_dur, inh_dur] = get_sniff_index(unit, md, params, sniff);
    
    ind = find(ind);
    if ~strcmp(params.r_type, 'rt') && ~strcmp(params.r_type, 'rwt')
        all = spikeAlignment(unit(md.unit_ref), snWF, ind, inh_dur, params);
        totalSniffs = size(unit(md.unit_ref).('rt'), 2);
        t_vec = size(all, 1);
        unit(md.unit_ref).(params.r_type) = sparse(zeros(t_vec, totalSniffs));
        unit(md.unit_ref).(params.r_type)(:,ind) = sparse(all);
    else
        all = unit(md.unit_ref).(params.r_type)(t_axis,ind);
    end
    
    
%     if nargin == 6 && strcmp(varargin{1}, 'by_inhal_dur')
%         [~, sn_ind] = sort(sn_dur);
%         [sp_pos, tr_num] = find(all(:, sn_ind));
%     else
%         [sp_pos, tr_num] = find(all);
%     end
    if nargin == 6 && strcmp(varargin{1}, 'by_sniff_dur')
        [sn_dur, sn_ind] = sort(sn_dur);
        [sp_pos, tr_num] = find(all(:, sn_ind));
    elseif nargin == 6 && strcmp(varargin{1}, 'by_inhal_dur')
        [inh_dur, inh_ind] = sort(inh_dur);
        [sp_pos, tr_num] = find(all(:, inh_ind));
    else
        [sp_pos, tr_num] = find(all);
    end
end

function spikes = spikeAlignment(un, sn, rind, dur, params)
% the function use spikes from rt matrix for each response and creates
% different spike alignment depend on option:
% options:
%  - warp_warp - warping separately of inhalation and exhalation
%  - time_warp - inhalation in time, exhalation warped
%  - volume - by the inhalation volume
%  - scaling - exhalation scaled proportionally to inhalation 
%
% output is a structure of rasters for each alignment
%
%                                       Roma               June 10 2016


[inh_dur, ind]  = sort(dur);
inh_ind = rind(ind);

%calculation of mean parameters
mean_inh  = round(mean(un.inh_dur));
mean_sn   = round(mean(un.sniff_dur));
mean_exh  = mean_sn - mean_inh;

% os - is for odor specific
os_mean_inh  = round(mean(un.inh_dur(rind)));


% rt   = zeros(max(ceil(un.sniff_dur(rind))), length(rind));
switch params.r_type
    case 'vol'
        spikes  = zeros(params.range(2), length(rind));
end
        
yind = 0;


% rectification of negative pressure part
sn.waveform(sn.waveform<0)=0;
% integral calculation for each sniff
maxVolume = zeros(1, length(inh_dur));
for k = 1:length(inh_dur)
    volume = cumsum(int32(sn.waveform(inh_ind(k), :)));
    maxVolume(k) = max(volume);
end

% search for time shift
mv = sort(maxVolume);
tshift = mv(end-3)/20;


        
% spike times realignment according to different models
for ir = inh_ind
    
    yind = yind+1;
    
    t1 = round(un.inh_dur(ir));
    t2 = ceil(un.sniff_dur(ir));
 
    if t2<=size(un.rt, 1)
        
        all_sp = find(full(un.rt(1:t2, ir))==1);
        switch params.r_type
            case 'vol'
                volume = cumsum(int32(sn.waveform(inh_ind(k), :)));
                sp = round(all_sp) + (os_mean_inh - t1) - find(volume>0.1*tshift, 1);
                spikes(sp(sp>0), yind) = 1;
%                 rvol((os_mean_inh+(t2-t1)+1):end, yind) = NaN;
        end

%         % 'time dilation'
%         rdil_spikes = round(all_sp/t1*mean_inh);
%         rdil(rdil_spikes, yind) = 1;
%         rdil((max(rdil_spikes)+1):end, yind) = NaN;
    end
end


end

function [psth, t_axis, mean_psth, std_psth] = binned_raster(unit, md, params, sniff, binsize)

    range       = params.range;
    r_type      = params.r_type;

    ind = get_sniff_index(unit, md, params, sniff);        

    rast = full(unit(md.unit_ref).(r_type)(range(1):range(2),ind)');

    ntr = size(rast, 1);
    np  = size(rast, 2);

    if rem(np,binsize)>0
        np = np - rem(np,binsize);
    end

    [y, x] = find(rast);
    psth = zeros(ntr, np/binsize);
    
    for i = 1:max(y)
        psth(i, :) = histcounts(x(y==i), 0:binsize:np)/binsize*1e3;
    end
    
    t_axis = binsize:binsize:np;
    mean_psth = mean(psth, 1);
    std_psth = std(psth, 1);
end

function [fast_sp_count, slow_sp_count, sn_dur] = spCountFastVsSlow(unit, md, params, sniff)
% builds raster matrix and if handle to axes provided, 
% it will build raster plot
%
% [all,t_axis, g] = dc.raster(unit, md(1), params, 1);
%
%                                           Roma       26 Feb 2016
    if length(md)>1;  error('one md at a time pls');  end

    range       = params.range;    
    t_axis      = range(1):range(2);
    
    [ind, sn_dur] = get_sniff_index(unit, md, params, sniff);
    all = unit(md.unit_ref).(params.r_type)(t_axis,ind);  
    
%     [~, sn_ind] = sort(sn_dur);
    [~, tr_num] = find(all(1:120, sn_dur(sn_dur<250)));
%     figure; plot(a, tr_num, '.')
    fast_sp_count(1:max(tr_num)) = accumarray(tr_num, 1);
    [~, tr_num]  = find(all(121:end, :));
%     figure; plot(a, tr_num, '.')
    slow_sp_count(1:max(tr_num)) = accumarray(tr_num, 1);
end

function [b, bint, lm] = regr(unit, md, params, sniff)
% performs regression analysis on given sniff cycle of activity vs. 
% sniff frequency.
% [all,t_axis, g] = dc.raster(unit, md(1), params, 1, g, 'by_freq');
%
%                                           Roma       26 Feb 2016

%     n = 100;
%     noise = randn(n, 1);
%     x = 10*rand(n, 1);
%     y = 2 + 3.5*x + noise;
%     X = [ones(size(x)) x];
%     plot(x, y, '.')
%     lsline

    [~, tr_num, sn_dur] = raster4regr(unit, md, params, sniff);
    sp_count = zeros(length(sn_dur), 1);
    sp_count(1:max(tr_num)) = accumarray(tr_num, 1);
    freq = 1000./sn_dur';
    plot(freq, sp_count, '.')
    lsline
   
    Fr = [ones(size(freq)) freq];
 
    [b,bint] = regress(sp_count, Fr);
    lm = fitlm(Fr,sp_count,'linear');
end

function [avg,t_axis,all,avg_resid,err,err_r] = psth(unit, md, params, sniff)
% returns the psth for a given md structure and sniff_num
% 
% params.part and params.n_split now manage what portion of the trials are
% returned in the psth
%
% if n_split is Inf then pick a random trial for each cell
% also compute the complement average excluding the selected trial and
% return as avg_resid
%
%                                           Roma        17 Feb 2016
    range       = params.range;
    wnd         = params.wnd;
    step        = params.step;
    r_type      = params.r_type;
    
    n_split     = params.n_split;
    part        = params.part;
    if part>n_split || part<0; error('part must be > 0 & < n_split'); end
    
    n_md        = length(md);
    nw          = ceil((diff(range)-wnd+2)/step);
    in_aver     = (1:wnd)'*ones(1,nw) + ones(wnd,1)*(0:step:diff(range)+1-wnd);
    t_axis      = (range(1)-1+wnd/2):step:(range(2)-wnd/2);    
    
    avg         = zeros(nw,n_md);
    avg_resid   = zeros(nw,n_md);
    
    err         = zeros(nw,n_md);
    err_r       = zeros(nw,n_md);    
    
    all = cell(1,n_md);
    for n = 1:n_md
        ind = get_sniff_index(unit, md(n), params, sniff);        
        ind = find(ind);
        
        if isinf(n_split)
            tmp = ind(randi(numel(ind),1,1)); % get a random trial
            ind_r = setdiff(ind,tmp); % get the complementary trials
            ind = tmp;
        else
            part = 1;
            ind = ind(part:n_split:end);
        end
        
        spk = unit(md(n).unit_ref).(r_type)(range(1):range(2),ind);
        all{n} = 1000*shiftdim(mean(reshape(full(spk(in_aver,:)),[wnd, nw, numel(ind)]),1),1);        
        avg(:,n) = mean(all{n},2);
        err(:,n) = std(all{n},[],2);
        
        if isinf(n_split)
            spk = unit(md(n).unit_ref).(r_type)(range(1):range(2),ind_r);
            tmp = 1000*shiftdim(mean(reshape(full(spk(in_aver,:)),[wnd, nw, numel(ind_r)]),1),1);
            avg_resid(:,n) = mean(tmp,2);
            err_r(:,n) = std(tmp,[],2);
        end
        
    end    
end

function [avg,t_axis,all,err, slow, slow_avg,...
          slow_err, fast, fast_avg, fast_err] = split_psth(unit, md, params, sniff)
      
% slow, slow_avg, slow_err, fast, fast_avg, fast_err

% returns the psth for a given md struture and sniff_num
% 
% params.part and params.n_split now manage what portion of the trials are
% returned in the psth
%
% if n_split is Inf then pick a random trial for each cell
% otherwise it splits trials by frequency.
%
% also compute the complement average excluding the selected trial and
% return as avg_resid
%
% 
%                                           Roma        17 Feb 2016
    range       = params.range;
    wnd         = params.wnd;
    step        = params.step;
    r_type      = params.r_type;

    n_split     = params.n_split;
    part        = 0.3;%params.part;
    if part>n_split || part<0; error('part must be > 0 & < n_split'); end

    nw          = ceil((diff(range)-wnd+2)/step);
    in_aver     = (1:wnd)'*ones(1,nw) + ones(wnd,1)*(0:step:diff(range)+1-wnd);
    t_axis      = (range(1)-1+wnd/2):step:(range(2)-wnd/2);

    [ind, sn_dur] = get_sniff_index(unit, md, params, sniff);
    [~, sn_ind] = sort(sn_dur);
    ind = find(ind);
    ind = ind(sn_ind); 

    ind = ind(1:n_split:end);

    splitInd = ceil(part*length(ind));
    indSlow = ind((splitInd+1):end);
    indFast = ind(1:splitInd);

    spk = unit(md.unit_ref).(r_type)(range(1):range(2),indFast);
    fast = 1000*shiftdim(mean(reshape(full(spk(in_aver,:)),[wnd, nw, numel(indFast)]),1),1);
    fast_avg = mean(fast,2);
    fast_err = std(fast,[],2)./size(fast,2);
    
    spk = unit(md.unit_ref).(r_type)(range(1):range(2),indSlow);
    slow = 1000*shiftdim(mean(reshape(full(spk(in_aver,:)),[wnd, nw, numel(indSlow)]),1),1);
    slow_avg = mean(slow,2);
    slow_err = std(slow,[],2);
    
    spk = unit(md.unit_ref).(r_type)(range(1):range(2),ind);
    all = 1000*shiftdim(mean(reshape(full(spk(in_aver,:)),[wnd, nw, numel(ind)]),1),1);
    avg = mean(all,2);
    err = std(all,[],2);
    
    kernel      = ones(step,1)/(step);
    slow_avg    = conv(slow_avg,kernel,'same');
    fast_avg    = conv(fast_avg,kernel,'same');

end

function [m_phase, r_phase] = get_phase(unit,md,params,sniff,plotit)
% Compute circular phase of response assuming 320 msec sniff cycle
%
%                                           Yevgeniy    01/09/12
    if nargin<5; plotit = 0; end
    m_sniff_dur = 500;
    
    m_phase = zeros(1,numel(md));
    r_phase = zeros(1,numel(md));
    
    for n = 1:numel(md)
        [r,t] = psth(unit,md(n),params,sniff);
        %r = r-mean(r);
                
        angle = (t/m_sniff_dur*2*pi)';      % the phase within the sniff
                
        c = [sin(angle) cos(angle)];        % convert to cartesian       
        c = mean(c.*[r r])./mean([r r]);    % compute weighted average
        m_phase(n) = atan2(c(1),c(2));      % mean phase of resultant
        r_phase(n) = sqrt(sum(c.^2));       % amplitude of mean vector
        
        if plotit                           % plot result
            figure;
            polar(angle,r,'k');
            hold on;
            polar(m_phase(n)*[1 1],r_phase(n)*[0 1],'r')
        end
    end
end

function doublets = get_doublets(md)
% return pairs of indices in md corresponding to two concentrations of the
% same odor presented to the same cell
% rows are each doublet
%
%                                           Yevgeniy        01.02.2012
    count = 0;
    doublets = [0 0];
    
    mice = unique([md.mouse]);
    for m = mice % for each mouse
        ind_m = [md.mouse] == m;
        sesss = unique([md(ind_m).sess]);
        for s = sesss % for each session
            ind_s = [md.sess] == s;
            units = unique([md(ind_s & ind_m).cell_num]);
            for u = units % for each cell of that session
                ind_u = [md.cell_num] == u;
                valves = unique([md(ind_u & ind_s & ind_m).valve]);
                for v = valves % for each valve
                    ind_v = [md.valve] == v;
                    concs = unique([md(ind_u & ind_s & ind_m & ind_v).conc]);
                    for c = 1:numel(concs)-1
                        ind_c1 = [md.conc] == concs(c);
                        ind_c2 = [md.conc] == concs(c+1);
                        count = count + 1;
                        
                        try
                            doublets(count,1) = find(ind_u & ind_s & ind_m & ind_v & ind_c1);
                            doublets(count,2) = find(ind_u & ind_s & ind_m & ind_v & ind_c2);
                        catch
                            keyboard
                        end
                    end
                end
            end
        end
    end
end

function triplets = get_triplets(md)
% return triplets of indices in md corresponding to three concentrations of the
% same odor presented to the same cell
% rows are each triplet
%
%                                           Yevgeniy        01.02.2012
    count = 0;
    triplets = [0 0 0];
    
    mice = unique([md.mouse]);
    for m = mice % for each mouse
        ind_m = [md.mouse] == m;
        sesss = unique([md(ind_m).sess]);
        for s = sesss % for each session
            ind_s = [md.sess] == s;
            units = unique([md(ind_s & ind_m).cell_num]);
            for u = units % for each cell of that session
                ind_u = [md.cell_num] == u;
                valves = unique([md(ind_u & ind_s & ind_m).valve]);
                for v = valves % for each valve
                    ind_v = [md.valve] == v;
                    concs = unique([md(ind_u & ind_s & ind_m & ind_v).conc]);
                    for c = 1:numel(concs)-2
                        ind_c1 = [md.conc] == concs(c);
                        ind_c2 = [md.conc] == concs(c+1);
                        ind_c3 = [md.conc] == concs(c+2);
                        count = count + 1;
                        
                        try
                            triplets(count,1) = find(ind_u & ind_s & ind_m & ind_v & ind_c1);
                            triplets(count,2) = find(ind_u & ind_s & ind_m & ind_v & ind_c2);
                            triplets(count,3) = find(ind_u & ind_s & ind_m & ind_v & ind_c3);
                        catch
                            keyboard
                        end
                    end
                end
            end
        end
    end
end

function [pop,pop_r,pop_err,pop_r_err] = pop_vector(unit,md,params,sniff,flag)
% return population vector for a given sniff
% flag = 0 : leave mean in
% flag = 1 : take mean out
% flag = 2 : keep only mean
%                                           Yevgeniy        01.02.2012    
    nomean = 0;
    if nargin>4
        nomean = flag;
    end
        
    if nomean==2
        params.wnd = diff(params.range)+1;
        params.step = 1;
    end
    
    n_md        = length(md);
    nw          = ceil((diff(params.range)-params.wnd+2)/params.step);
    
    pop         = zeros(n_md*nw,numel(sniff));
    pop_r       = zeros(n_md*nw,numel(sniff));

    pop_err     = zeros(n_md*nw,numel(sniff));
    pop_r_err   = zeros(n_md*nw,numel(sniff));

    for n = 1:numel(sniff)
        [tmp,~,~,tmp_r,err,err_r] = psth(unit,md,params,sniff(n));
        if nomean==1
            tmp =     tmp - repmat(mean(tmp  ,1),[nw 1]);
            tmp_r = tmp_r - repmat(mean(tmp_r,1),[nw 1]);
        end
        pop(:,n)   = reshape(tmp,  [n_md*nw,1]);      
        pop_r(:,n) = reshape(tmp_r,[n_md*nw,1]);     
        
        pop_err(:,n)   = reshape(err,  [n_md*nw,1]);      
        pop_r_err(:,n) = reshape(err_r,[n_md*nw,1]);
    end
end

function [cxy,lags,dt,mc] = md_xcorr(unit,md,params,sniffs)
% retrieve the cross correlation function for average PSTH between
% two md structs
%
% md should specify two cell-odor pairs to compare
% sniffs specifies which sniffs to compare for each cell odor pair
%
%                                           Yevgeniy        01.02.2012
%
    if numel(md)~=2; error('Need exactly 2 md structs'); end;       
    nw          = ceil((diff(params.range)-params.wnd+2)/params.step/2);
    
    v1 = psth(unit,md(1),params,sniffs(1));
    v2 = psth(unit,md(2),params,sniffs(2));
    [cxy,lags] = xcorr(v1-mean(v1),v2-mean(v2),nw,'coeff');
    [mc,pos] = max(cxy);
    lags = lags*params.step;
    dt = lags(pos);
end

function [cxy,lags,dt,mc] = md_xcorr_boot(unit,md,params,sniffs)
% retrieve the cross correlation function PSTH between
% two md structs with bootstrap
%
% md should specify two cell-odor pairs to compare
% sniffs specifies which sniffs to compare for each cell odor pair
%
%                                           Yevgeniy        01.09.2012
%
    if numel(md)~=2; error('Need exactly 2 md structs'); end       
    nw          = ceil((diff(params.range)-params.wnd+2)/params.step/2);
    
    [~,~,all1] = psth(unit,md(1),params,sniffs(1));
    all1 = all1{1};
    v2 = psth(unit,md(2),params,sniffs(2));
    v2 = v2-mean(v2);
    
    for n = 1:size(all1,2)
        % for each trial of the first md, compare to 2nd md mean
        %ind = n;
        ind = randsample(size(all1,2),floor(size(all1,2)/2));
        
        v1 = mean(all1(:,ind),2);
        v1 = v1-mean(v1);

        [cxy(:,n),lags] = xcorr(v1,v2,nw,'coeff');
        [mc(n),pos(n)] = max(cxy(:,n));
    end
    lags = lags*params.step;
    dt = lags(pos);
end

function class_success = cell_classify(Y, X, nruns)
%%%% single cell classification %%%%%%%%%%%%%
% As input it gets pair of concentrations and output is a vector of 
% cumulative discrimination success
%                                           Roma        10 May 2016
%
    [ntr1, nbins] = size(Y);
%     nbins = 1;
    ntr2 = size(X, 1);
  
    P = zeros(nruns, nbins);
    for r = 1:nruns
        rand_ind = randi(ntr2); % pulling out random trial
        trial = X(rand_ind, :);
        
%         test1 = Y(randi(ntr1), :); % test1 consist of 70% of trials from Y
%         test2 = X(randi(ntr2), :); % test2 consist of 70% of trials from X
        
        test1 = mean(Y(randi(ntr1, round(0.7*ntr1), 1), :), 1);
        test2 = mean(X(randi(ntr2, round(0.7*ntr2), 1), :), 1);
        
        corr = zeros(1, nbins);
        for ib = 1:nbins
            %calculating euclidean distance to the centers of two groups  
            dist2wrongClust  = sum((test1(:, 1:ib)-trial(1:ib)).^2);
            dist2corrClust = sum((test2(:, 1:ib)-trial(1:ib)).^2);
            if dist2corrClust==dist2wrongClust
                corr(ib) = randi(2,1,1);
            else
                [~, corr(ib)] = min([dist2wrongClust, dist2corrClust]);
            end
            
        end
        corr = corr-1; % after substraction 1 is for correct group and 0 for incorrect
        P(r, :) = corr;      
    end
    class_success = mean(P, 1);
    
%     if max(class_success)>0.5
%         class_success = 1 - class_success;
%     end
end

function [p_match,D] = pop_classify(unit,md,params,doublets,cL,cR,flag)
% classify random samples from responses to different odors/sniffs
%
%                                           Yevgeniy        01.12.2012
%

    N_SAMPLES   = 100;    

    train = [];
    train_err = [];
    
    n_cR = size(cR,1);
    n_cL = size(cL,1);
    
    D = zeros(N_SAMPLES,n_cR,n_cL);
    
    % Construct training set from cR
    for r = 1:n_cR
        [train(:,r),~,train_err(:,r),~] = pop_vector(unit,md(doublets(:,cR(r,1))),params,cR(r,2),flag);
    end
    
    params.n_split  = Inf;
    
    for n = 1:N_SAMPLES
        % for each response type to be classified                
        for l = 1:n_cL
            train_use = train;
            train_err_use = train_err;

            % Draw a sample from cL
            [sample,rest,~,rest_err] = pop_vector(unit,md(doublets(:,cL(l,1))),params,cL(l,2),flag);

            % replace any matching samples if needed
            pos = find(sum(repmat(cL(l,:),n_cR,1)==cR,2)==2,1,'first');
            if ~isempty(pos)
                train_use(:,pos) = rest;
                train_err_use(:,pos) = rest_err;
            end
            
            %err = repmat(sqrt(sum(train_err_use.^2,2)),1,n_cR);            
            %D(n,:,l) = sqrt(sum(((train_use-repmat(sample,1,n_cR))./err).^2));
            sample = sample - mean(sample);            
            train_use = train_use-repmat(mean(train_use),size(train_use,1),1);
            D(n,:,l) = sqrt(sum(((train_use-repmat(sample,1,n_cR))).^2));              
            
        end
    end
    
    for n = 1:n_cR
        [~,match] = min(D(:,:,n),[],2);
        p_match(:,n) = hist(match,1:n_cR);
    end
end

function [cxy,lags,n_lag,p_match_corr,p_match_euclid] = pop_xcorr_boot(unit,md,params,doublets,cL,cR,match_flag,flag)
% compute cross correlation for 10000 random "simulated" trials compared to
% average "template" vectors
% 
% cL - indicates the concentration/sniff pair for left hand input into xcorr
%      this is the set of responses from which a single trial is drawn
%
% cR - indicates the concentration/sniff pair for right hand input into xcorr
%    - this is the set of responses from which the basis is drawn
%
% n_lag - the number of a given lag observed in cxy
%
% p_match - probability of matching a given cL to a given cR
%
%                                           Yevgeniy        01.05.2012
%
    n_L             = size(cL,1);
    n_R             = size(cR,1);
    n_match         = sum(match_flag);
    N_BOOT          = 10000;
    MAX_LAG         = 100;
    lags            = -MAX_LAG:MAX_LAG;
    
    if flag == 2; 
        MAX_LAG     = 0; 
        lags        = 0;
    end;
    
    vecR = cell(1,n_R);
    for n = 1:n_R
        vecR{n} = pop_vector(unit,md(doublets(:,cR(n,1))),params,cR(n,2),flag);
    end
    
    params.n_split  = Inf;

    cxy             = zeros(2*MAX_LAG+1,N_BOOT,n_R,n_L);
    C               = zeros(N_BOOT,n_R,n_L);
    E               = zeros(N_BOOT,n_R,n_L);
    
    n_lag           = zeros(2*MAX_LAG+1,n_R,n_L);
    p_match_corr    = zeros(n_match,n_L);
    p_match_euclid  = zeros(n_match,n_L);
    
    for L = 1:n_L 
        for n = 1:N_BOOT            
            [vecLt,temp] = pop_vector(unit,md(doublets(:,cL(L,1))),params,cL(L,2),flag);       
            
            for R = 1:n_R
                vecRt = vecR{R};                                
                if sum(cL(L,:)==cR(R,:))==2; vecRt = temp; end;
                
                cxy(:,n,R,L) = xcorr(vecLt,vecRt,MAX_LAG,'coeff');
                C(n,R,L) = corr(vecLt,vecRt);
                E(n,R,L) = sum(abs(vecLt-vecRt));
            end
        end
        
        % Compute lags
        if flag == 2
            n_lag = [];
        else
            for R = 1:n_R
                [~,lag] = max(cxy(:,:,R,L),[],1);
                n_lag(:,R,L) = hist(lags(lag),lags);
            end
        end
        fprintf('Completed %d of %d...\n',L,n_L);
       
        % Compute match using correlation     
        [~,match] = max(C(:,match_flag,L),[],2);
        p_match_corr(:,L) = hist(match,1:n_match);
        p_match_corr(:,L) = p_match_corr(:,L)/N_BOOT;

        [~,match] = min(E(:,match_flag,L),[],2);
        p_match_euclid(:,L) = hist(match,1:n_match);
        p_match_euclid(:,L) = p_match_euclid(:,L)/N_BOOT;

    end
    
end

function out = has_response(unit, md, params, sniff)
% detect statistically significant deviations from background activity
%
%                                           Roma        01.06.2012
%

    params.wnd      = 10; % switch window to 30 msec (most informative)
    params.step     = 10; % high step size to decrease number of comparisons
    n_md            = length(md);
    nw              = ceil((diff(params.range)-params.wnd+2)/params.step);
    N_BOOT          = 20000;
    
    bckg_mean_min   = 10;
    amp_coef        = 0.5;
    out             = md;
            
    fig = figure;
    for n = 1:n_md
        % calculation of average psth of 3 cycles prior odor onset 
        if (n==1) || (md(n).unit_ref~=md(n-1).unit_ref)
            % fractional number for sniff means sniff prior odor onset
            [avg_0a1,~,all_01] = psth(unit, md(n), params, 0.1);
            all_01 = all_01{1};
            [avg_0a2,~,all_02] = psth(unit, md(n), params, 0.2);
            all_02 = all_02{1};
            [avg_0a3,~,all_03] = psth(unit, md(n), params, 0.3);
            all_03 = all_03{1};
            avg_0a = (avg_0a1 + avg_0a2 + avg_0a3)/3;
            all_0 = [all_01 all_02 all_03];
            %all_0 = all_0{1};
            n_sniff0 = size(all_0,2);            
        end
        
        [avg_1,~,all_1] = psth(unit, md(n), params, sniff);
        all_1 = all_1{1};
        n_sniff1 = size(all_1,2);        

        out(n).resp_latency = 0;
        out(n).has_response = 0;
        out(n).initial_sign = 0;
                
        % 1st criterion: mean(background) > 10
        crit1 = mean(avg_0a) > bckg_mean_min;
        
        % 2nd criterion: change in firing rate > mean(background)      
        ind = find(abs(avg_1-mean(avg_0a))/mean(avg_0a) > amp_coef, 1, 'first');
        crit2 = ~isempty(ind) & crit1;
        
        if crit2
            
            % question is how likely are we to get this avg1 firing rate by
            % chance from all_0
            
            fprintf('Cell %s Boot N=%d (pop = %d)...',md(n).ID, n_sniff1, n_sniff0);
            
            % sample n_sniff1 trials from all_0 N_BOOT times and examine
            % distribution
            avg_0 = zeros(nw,N_BOOT);
            for b = 1:N_BOOT
                %ind = randperm(n_sniff0);
                ind = randsample(n_sniff0,n_sniff1);
                avg_0(:,b) = mean(all_0(:,ind(1:n_sniff1)),2);
            end
            
            % find points outside 99% confidence interval corrected for number
            % of comparisons
            fprintf('...p_adj = %1.4f (%d)...',0.01/nw,floor(N_BOOT*0.005/nw));
            
            avg_0 = sort(avg_0,2);
            avg_0_lo = avg_0(:,floor(N_BOOT*0.005/nw));
            avg_0_hi = avg_0(:,ceil(N_BOOT*(1-0.005/nw)));
            
            sig_lo = avg_1<avg_0_lo;
            sig_hi = avg_1>avg_0_hi;
            
            p_lo = find(sig_lo,1,'first');
            p_hi = find(sig_hi,1,'first');
            
            figure(fig);
            plot(avg_0_lo);
            hold on;
            plot(avg_0_hi)
            plot(avg_0a);
            plot(avg_1,'r');
            hold off;
            if ~isempty([p_lo p_hi])
                
                fprintf(' ... RESPONSIVE!\n');
                drawnow;
                %pause;
                
                out(n).has_response = 1;
                if ~isempty(p_lo) && (isempty(p_hi) || (p_lo<p_hi))
                    out(n).initial_sign = -1;
                    out(n).resp_latency = p_lo;
                elseif ~isempty(p_hi) && (isempty(p_lo) || p_hi<p_lo)
                    out(n).initial_sign = 1;
                    out(n).resp_latency = p_hi;
                end
            else
                fprintf(' ... no\n');
            end
        else
            fprintf('Cell %s CRIT1 FAIL\n',md(n).ID);
        end
    end
end

function out = sharp_responses(unit, md, sniff_num, params)
% detect sharp events
%                                           Roma            12.28.2011
% function detect sharp responses for md and given sniff_num, 
% if sniff_num is not specified sniff_num = 1;
% if params are not specified params = default_params();
%
% the criteria for sharp response detection are the following:
% 1) mean(bckg) > bckg_mean_min          parameter: bckg_mean_min = 10
% 2) there is exist time point for which rate(t) > amp_coef*bckg(t)
%                                        parameter: ampl_coef = 3
% 3) for such temporal point preceeding firing rate increase is sharp: 
%    rate(t) > incr_coef*rate(t-delta)   parameter: incr_coef = 2, 
%                                                   delta = 3*wnd/step
% additional parameters:
% isi_coef         ratio between inverted rate peak amplitude and minimal
%                  isi for sharp event detection
% se_pre, se_post: interval around peak of firing rate where sharp even is
%                  detected
% to_print:        to print the raster plots and psth
% to_pause:        to pause at the figures with sharp events
%
% the function returns the structure which fields are identical to md
% strucutre with addition of the fields se (sharp_response):
% se_exits: 0, 1 
% se_lat:   sharp response latency, the mean latency of the first spike
%           preceeding small ISI spieks
% se_rel:   sharp response reliability (the proportion of trials for which
%           sharp increase of firing rate exist
% se_jitter: sharp response jitter, std of the latency of the first spikes
%           preceeding spike swith small ISI
        
    if ~exist('params', 'var')
        params = default_params();
    end
    if ~exist('sniff_num', 'var')
        sniff_num = 1;
    end
    
    range       = params.range;
    wnd         = params.wnd;
    step        = params.step;

    % sharp response parameteres
    amp_coef      = 3;
    bckg_mean_min = 10;
    incr_shift    = ceil(3*wnd/step);
    incr_coef     = 2;
    isi_coef      = 1.5;
    se_pre        = -2;
    se_post       = 2;
    to_print      = 1;
    to_pause      = 0;

    out    = md;
    ku_cur = 0;
    
    for id = 1:length(md)  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ku = md(id).unit_ref;
        
        if ku ~= ku_cur
            ku_cur = ku;
            [bckg, t_axis] = psth(unit, md(id), params, 0);
        end
            
        se_rel   = 0;
        se_lat   = 0;
        se_jit   = 0;
        
        % psth estimation
        rate  = psth(unit, md(id), params, sniff_num);
        rr    = raster(unit, md(id), params, sniff_num);
        n_tr  = size(rr,2);

        if to_print    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            figure(1), clf
            gs1 = subplot(2,1,1);
            [sp, tr] = find(rr);
            plot(sp, tr, '.k')
            xlim(range)
            title(sprintf('un: %d   od: %d    conc: %4.2f', ku, md(id).valve, md(id).conc))
     
            gs2 = subplot(2,1,2);
            plot(t_axis, bckg, 'k'), hold on
            plot(t_axis, rate, 'b'), hold off
        end    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        fprintf('unit: %d   odor: %d   conc: %4.2f  ', ku, md(id).valve, md(id).conc);
                
        % 1st criterium: mean(bckg) > 10
        se_exist = mean(bckg) > bckg_mean_min;
        fprintf(' 1: %d   ', se_exist)
        
        % 2nd criterium: rate(t) > 3*bckg(t)
        if se_exist
            ind = find(rate > amp_coef*mean(bckg), 1, 'first');
            se_exist = ~isempty(ind);
            fprintf(' 2: %d   ', se_exist)
        end
        
        % 3rd criterium:  psth(t) > 3*rate(t-delta)
        if se_exist
            ind1 = max([1,ind-incr_shift]);
            se_exist = rate(ind) > incr_coef*rate(ind1);
            fprintf(' 3: %d   ', se_exist)
        end

        if se_exist && id~=361 && id~=362    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % looking for maximum in the vicinity of ind
            imax  = find(diff(rate(ind:end))<0, 1, 'first')-1+ind;
            % peak postions and isi
            peak_amp = rate(imax);
            peak_pos = t_axis(imax); 
%             if ku==41 && md(id).valve == 6
%                 keyboard
%             end
            isi_crit = ceil(1/peak_amp*1000*isi_coef);
               
            % boundaries for sharp even search 
            t1 = peak_pos + se_pre*isi_crit;
            t2 = peak_pos + se_post*isi_crit;
                                  
            if to_print    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                subplot(gs2),   hold on
                plot(peak_pos, peak_amp, 'or')
                plot([1;1]*[t1,t2], [0;peak_amp]*[1,1], 'c')
                hold off
                         
                subplot(gs1),   hold on
                plot([1;1]*[t1,t2], [0;n_tr+1]*[1,1], 'c')
                hold off
            end    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
            tr_se_exist = zeros(1, n_tr);
            tr_se_time  = zeros(1, n_tr);
            
            for it = 1:n_tr    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                % trial spikes in the search interval
                sp_int = sp((sp > t1)&(sp < t2)&(tr==it));
                % interspike interval
                isi = diff(sp_int);
                % find isi < isi crit
                in_isi = find(isi <= isi_crit);

                if ~isempty(in_isi)
                    tr_se_exist(it) = 1;
                    tr_se_time(it)  = sp_int(in_isi(1));
                            
                    if to_print
                        subplot(gs1), hold on
                        plot(sp_int(in_isi), it*ones(1,numel(in_isi)), '.r')
                        plot(sp_int(in_isi(1)), it, 'og')
                    end
                end
            end    % it ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            if any(tr_se_exist)>0
                se_rel = mean(tr_se_exist);
                se_lat = mean(tr_se_time(tr_se_exist > 0));
                se_jit = std(tr_se_time(tr_se_exist > 0));
                      
                if to_print
                    subplot(gs2), hold on
                    plot(se_lat*[1,1], [0, peak_amp], 'r', 'LineWidth',2)
                    plot(se_lat+se_jit*[-1,1], peak_amp*[1,1], 'r', 'LineWidth',2)
                end
            else 
                se_exist = 0;
            end
            
            fprintf('   lat: %4.1f   prec: %3.2f   reliab:%3.2f', se_lat, se_jit, se_rel)
            if to_pause
                pause
            end
        end    % se_exist ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        
        out(id).se_exist  = se_exist;
        out(id).se_lat    = se_lat;
        out(id).se_rel    = se_rel;
        out(id).se_jitter = se_jit;
        
        fprintf('\n')
    end    % id ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end

function [unit, sn] = all_units()
% builds unit structure for all sessions
%                                       Dima                12.29.2011
    SESS = {[ 57,26],   [ 57,32],...
            [ 77, 1],   [ 75, 2],...
            [553, 2],... % sesison [075, 1] is missing sniff
            [576, 1],   [576, 2],...
            [581, 1],   [581, 2],...
            [583, 2],   [ 77, 4],...
            [104, 1],   [105, 3],...
            ... % following sessions were aquuires with Evgeniy
            [110, 5],   [110, 6], ...
            [111, 3],   [111, 6], ...
            [113, 2],   [113, 3], ...
            [114, 2],   [114, 3], ...
            [114, 4],   [114, 5]};
        
    for is = 1:length(SESS)
        fprintf('Creating unit for m_%03d_s%02d...\n',SESS{is});
        if is == 1
            unit = create_unit(SESS{1}(1), SESS{1}(2));
            sn = create_sniff(SESS{1}(1), SESS{1}(2));
        else
            unit = [unit, create_unit(SESS{is}(1), SESS{is}(2))];
            sn = [sn, create_sniff(SESS{is}(1), SESS{is}(2))];
        end
    end
            

end

function un = create_unit(mouse, sess)
% builds unit structure for a given mouse and session
%                                   Dima                    12.29.2011
%
% required files:   sess_mmm_ss.mat, spikes_mmm_ss.mat, sniff_mmm_ss.mat
% output strucutre fileds
%  mouse           - mouse number
%  sess            - session number
%  cell_num        - unit number in spikes_mmm_ss.mat strucutre
%  sorted          - sorted field in spikes_mm_ss.mat
%  odorID          - odorID from sess_info
%  t0              - time of the sniff onset             [1, n_sn] 
%  inh_dur         - inhalation duration                 [1, n_sn] 
%  sniff_dur       - sniff duration                      [1, n_sn] 
%  valve           - valve number                        [1, n_sn] 
%  conc            - concentration                       [1, n_sn] 
%  sn_odorSniffNum - odor sniff number                   [1, n_sn] 
%  rt              - spike rasters for each sniff sparse [nt,n_sn]; 
% 
    params = default_params();
    nt = params.nt;
    
    % fn - structure of file names for a given mouse and session
    fn = file_names(mouse, sess);
        
    % loading sniff
    disp(fn.sniff_waveforms);
    q  = load(fn.sniff_waveforms);
    sniff = q.sniff;
        
    % loading spikes
    disp(fn.db_spikes)
    q = load(fn.db_spikes);
    unit = q.unit;
        
    % load_sess
    disp(fn.db_sess)
    q      = load(fn.db_sess);
    trial  = q.trial;
    odorID = q.info.odorID;

    
    %----------------------------------------------------------------------
    % sniff data processing
    fprintf('sniff processing    number of sniffs: %d \n', length(sniff))
    n_sn = length(sniff);
    n_tr = length(trial);
    t_zer     = [sniff.t_zer];
    t_zer_fit = reshape([sniff.t_zer_fit],4, n_sn) - ones(4,1)*t_zer(1,:);
    tr_t0     = [trial.start];
    tr_od     = [trial.odorTime];
    
    % make compatible with previous versions
    if size(t_zer,1)==2
        warning('Only two entries in t_zer!','ac_tools');
        t_zer(3,:) = t_zer(2,:);
        t_zer(2,:) = t_zer_fit(2,:);
    end
 
    % time since the beginning of the session
    % correction for 32 msec
    sn_t0     = [sniff.t0] + t_zer(1,:) - params.correction;        
    
    % sniff duration
    sn_dur    = t_zer(3,:) - t_zer(1,:);
    
    % inhalation duration
    sn_inh    = t_zer_fit(2,:);
    % correction for wrong parabola fit:
    ind         = (imag(sn_inh) ~=0)|(sn_inh <=0)|(sn_inh >= sn_dur);
    sn_inh(ind) = t_zer(2,ind) - t_zer(1,ind);
    
    
    % building sn_trial array
    matr = ones(n_tr,1)*sn_t0 > tr_t0'*ones(1,n_sn);
    sn_trial = max(matr.*((1:n_tr)'*ones(1,n_sn)));
 
    % valve and conc sniff array
    matr = (ones(n_tr,1)*sn_t0 > (tr_t0 + tr_od(1,:))'*ones(1,n_sn))&(ones(n_tr,1)*sn_t0 < (tr_t0 + tr_od(2,:))'*ones(1,n_sn));
    sn_odorValve = sum(matr.*([trial.odorValve]'*ones(1,n_sn)),1);
    sn_odorConc  = sum(matr.*([trial.odorConc]'*ones(1,n_sn)),1);

    % odor sniff number array: 0 - bckg, 1,2,3 - sequential sniffs during 
    % odor presenation, -1,-2,-3 sniffs after odor presentation
    cur_num  = 0;
    sn_odorSniffNum = zeros(1, n_sn);
    for ks = 1:n_sn
        if sn_odorValve(ks) > 0
            cur_num = cur_num + 1;
        else
            switch sign(cur_num)
                case -1,    cur_num = cur_num - 1;
                case 0,     cur_num = 0;
                case +1,    cur_num = -1;
            end
            if cur_num < -3
                cur_num = 0;    
            end
        end
        sn_odorSniffNum(ks)  = cur_num;        
    end
    
    % unit strucutre
    ku      = 0;
    nt_max  = ceil(sn_t0(end))+nt;
    ind_sp  = (1:nt)'*ones(1,n_sn) + ones(nt,1)*ceil(sn_t0)-1;
    un      = struct();

    for iu = find([unit.sorted]>0)
        ku = ku + 1;
        un(ku).sess_ID      = sprintf('%03d_%02d', mouse, sess);    % session string ID
        un(ku).mouse        = mouse;                % mouse
        un(ku).sess         = sess;                 % session
        un(ku).cell_num     = iu;                   % unit number in spikes_mmm_ss.mat strucutre
        un(ku).sorted       = unit(iu).sorted;      % sorted field in spikes str
        un(ku).odorID       = odorID;               % odorID from sess_info
        un(ku).t0           = sn_t0;                % time of the sniff onset
        un(ku).inh_dur      = sn_inh;               % inhalation duration
        un(ku).sniff_dur    = sn_dur;               % sniff duration
        un(ku).valve        = sn_odorValve;         % valve number
        un(ku).conc         = sn_odorConc;          % concentration
        un(ku).odorSniffNum = sn_odorSniffNum;      % odor sniff number 
 
        % spike reformat
        rr = zeros(1, nt_max);
        rr(ceil(unit(iu).spikeTimes)) = 1;
        un(ku).rt = sparse(rr(ind_sp));
        
        fprintf('unit: %d \n', ku)
            
    end

end

function sn = create_sniff(mouse, sess)
% builds sniff structure for a given mouse and session
%                                   Roma                    11.27.2016

    
    % fn - structure of file names for a given mouse and session
    fn = file_names(mouse, sess);
        
    % loading sniff
    disp(fn.sniff_waveforms);
    q  = load(fn.sniff_waveforms);
    sniff = q.sniff;
        
        
    
    %----------------------------------------------------------------------
    % sniff data processing
    fprintf('sniff processing    number of sniffs: %d \n', length(sniff))
    n_sn = length(sniff);
    t_zer     = [sniff.t_zer];
    t_zer_fit = reshape([sniff.t_zer_fit],4, n_sn) - ones(4,1)*t_zer(1,:);

    
    % make compatible with previous versions
    if size(t_zer,1)==2
        warning('Only two entries in t_zer!','sd_tools');
        t_zer(3,:) = t_zer(2,:);
        t_zer(2,:) = t_zer_fit(2,:);
    end   
    
    % sniff duration
    sn_dur    = t_zer(3,:) - t_zer(1,:);
    sn.mouse = mouse;
    sn.sess = sess;
    sn.waveform = zeros(n_sn, 300, 'int16');
    sn_dur(sn_dur>=300)=299;
    for isn = 1:n_sn
        sn.waveform(isn,1:(sn_dur(isn)+1)) = int16(sniff(isn).waveform(t_zer(1,isn):t_zer(1,isn)+sn_dur(isn)));    
    end
end

function warp_par = mean_warp_par(unit, params)
% estimation of average inhalation and sniff duration for all individual
% sessions for sniffs in sniff_range
%                                   Dima                12.29.2011
 
    if ~exist('params', 'var')
        params = default_params();
    end

    sniff_range = params.sniff_range;
    
    sess_ID = [unit.mouse]*1000 + [unit.sess];
    sess = unique(sess_ID);
    ns   = length(sess);
    
    inh_dur   = zeros(1,ns);
    sniff_dur = zeros(1,ns);
    
    for is = 1:length(sess)
        ku  = find(sess_ID == sess(is), 1, 'first');
        ind = (unit(ku).sniff_dur >= sniff_range(1))&(unit(ku).sniff_dur <= sniff_range(2));
        inh_dur(is)   = mean(unit(ku).inh_dur(ind));
        sniff_dur(is) = mean(unit(ku).sniff_dur(ind));
    end
    
    warp_par = [mean(inh_dur), mean(sniff_dur)];
    
end

function un = warping(un, warp_par)
% the function use spikes from rt matrix for each unit and creates rw
% sparse matrix for each units.
% rw - is raster of warped spikes.
% spikes are warped to mean inhalation duration (warp_par(1)) and mean 
% sniff length (warp_par(2)),   2-interval warping
%                                       Dima                12.29.2011
%
% Fixed bug replacing:
%       rw(ceil(warp{ks}(rt_sn))) = 1;
% with:
%       rw(ceil(warp{ks}(rt_sn)),ks) = 1;
%                                       Yevgeniy            01.05.2012
% 
% rwt  - is a raster of half (inhalation) warped and rest 
%        (pause + exhal) is in time
%                                       Roma                04.23.2016

    params = default_params();
    nt = params.nt;
    
    T1 = warp_par(1);
    T2 = warp_par(2);
        
    sess_id   = [un.mouse]*1000 + [un.sess];
    cur_sess  = 0;
    warp      = cell(1,5e4);
    warp_time = cell(1,5e4);
    ind_sn    = cell(1,5e4);

    for iu = 1:length(un)  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if sess_id(iu) ~= cur_sess
            cur_sess = sess_id(iu);
            new_sess = true;
            n_sn = length(un(iu).t0);
            nt   = size(un(iu).rt, 1);
            disp('------------------------------------------')
        else
            new_sess = false;
        end
        
        rw   = zeros(T2, n_sn);
        rwt  = zeros(nt, n_sn);
        
        for ks = 1:n_sn    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if new_sess
                t1 = round(un(iu).inh_dur(ks));
                t2 = ceil(un(iu).sniff_dur(ks));

                warp{ks} = zeros(1,t2);
                warp{ks}(1:t1)    = (1:t1)/t1*T1;
                warp{ks}(t1+1:t2) = (1:t2-t1)/(t2-t1)*(T2-T1) + T1; 
                
                warp_time{ks} = zeros(1,nt);
                warp_time{ks}(1:t1)    = (1:t1)/t1*T1;
                warp_time{ks}(t1+1:nt) = (1:nt-t1) + T1; 
                
                ind_sn{ks} = 1:min([nt, un(iu).sniff_dur(ks)]);
            end
            % time index of spikes for unit iu and sniff ks
            rt_sn = un(iu).rt(ind_sn{ks},ks)==1;
            % warping
            rwt(ceil(warp_time{ks}(rt_sn)),ks) = 1;
            rw(ceil(warp{ks}(rt_sn)),ks) = 1;
        end    % ks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        un(iu).rw  = sparse(rw);
        un(iu).rwt = sparse(rwt);
        
        fprintf('unit: %d   sniffs: %d \n', iu, n_sn)
        
    end % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end

function [realTrace, timeModel, warpModel, warpTimeModel, volumeModel, timeDilationModel,...
          rt, rw, rwt, rvol, rdil] = spikeAlignmentModel(un, md, sn, bin)

% 1. Combine spikes during 1st sniff cycle into a vector and create a
% histogram.
% 2. For each model create PSTH and then using this PSTH create average
% spike vector.
%
%  - warp_warp - warping separately of inhalation and exhalation
%  - time_warp - inhalation in time, exhalation warped
%  - odor_dep_ww - odor specific warping of two pars
%  - odor_dep_tw - odor specific
%  - scaling - exhalation scaled proportionally to inhalation 
%
%                                   Roma               July 6 2016

my = my_tools();
maxPerc = 50;

% looking for the first sniff cycles for a given odor and given concentration
rind = find (un.valve == md.valve & un.conc == md.conc & un.odorSniffNum == 1);

[inh_dur, ind]  = sort(un.inh_dur(rind));
inh_ind = rind(ind);

%calculation of mean parameters
mean_inh  = round(mean(un.inh_dur));
mean_sn   = round(mean(un.sniff_dur));
mean_exh  = mean_sn - mean_inh;

% os - is for odor specific
os_mean_inh  = round(mean(un.inh_dur(rind)));
os_mean_sn  = round(mean(un.sniff_dur(rind)));
os_mean_exh = os_mean_sn - os_mean_inh;

% Estimation of the volume percentage that best aligns data
[volumes, r2vec] = bestL_Estimation(un, md, sn, bin, maxPerc);
[~, bestVolInd] = max(r2vec);
sprintf('best L - %02d%',volumes(bestVolInd))

% calculation of time shift for each cycle for specific percentage 
% each step is 5%
correction = timeShiftByVolume(sn, inh_ind, inh_dur, volumes(bestVolInd));
mean_corr = round(mean(correction));
fprintf('mean_corr - %3d\n', mean_corr)

% spike times realignment according to different models
[rt, rw, rwt, rvol, rdil]  = spikeRealignment(un, rind, correction);

%%%%%% time model %%%%%%
[rt_psth, ~]     = my.nan_spike_count(rt', bin); % in previous version I used nanpsth
tModel = reshape(repmat(rt_psth, bin, 1),[], 1);

%%%%%% warp model %%%%%%
[rw_psth, ~]     = my.nan_spike_count(rw', bin);
wModel = reshape(repmat(rw_psth, bin, 1), [], 1);

%%%%%% warp-time model %%%%%%
[rwt_psth, ~]    = my.nan_spike_count(rwt', bin);
wtModel = reshape(repmat(rwt_psth, bin, 1), [], 1);

%%%%%% volume model %%%%%%
[rvol_psth, ~]    = my.nan_spike_count(rvol', bin);
volModel = reshape(repmat(rvol_psth, bin, 1), [], 1);

%%%%%% time-dilation model %%%%%%
[rdil_psth, ~]  = my.nan_spike_count(rdil', bin);
tdModel = reshape(repmat(rdil_psth, bin, 1), [], 1);

realTraceDuration   = sum(un.sniff_dur(rind));
sniffTrace          = zeros(realTraceDuration, 1);
realTrace           = zeros(realTraceDuration, 1);
timeModel           = zeros(realTraceDuration, 1);
warpModel           = zeros(realTraceDuration, 1);
warpTimeModel       = zeros(realTraceDuration, 1);
volumeModel         = zeros(realTraceDuration, 1);
timeDilationModel   = zeros(realTraceDuration, 1);

vec = 1:length(inh_ind);
sn_ind = inh_ind;
t_ind = 0;
i = 0;
for ir = sn_ind
    i = i+1;
    t1 = ceil(un.inh_dur(ir));
    t2 = ceil(un.sniff_dur(ir));
    t2 = min(t2-rem(t2, bin), 900);

    %%%%%% sniff trace %%%%%%
    sniffTrace(t_ind+(1:300)) = sn.waveform(ir, :);
    
    %%%%%% real trace %%%%%%
    realTrace(t_ind+1:(t_ind+t2)) = full(un.rt(1:t2, ir))';
    
    %%%%%% time model %%%%%%
    timeModel(t_ind+1:(t_ind+t2)) = tModel(1:t2);
    
    %%%%%% warp model %%%%%%
%     interp_wModel = interp1(1:length(wModel), wModel, length(wModel)/t2:length(wModel)/t2:length(wModel), 'spline');
%     warpModel(t_ind+1:(t_ind+length(interp_wModel))) = interp_wModel;
    
    %%%%%% warp-time model %%%%%%
    % check if condition for 't_end' is correct
%     inhal_wtModel = interp1(1:os_mean_inh, wtModel(1:os_mean_inh), os_mean_inh/t1:os_mean_inh/t1:os_mean_inh , 'spline');
%     t_end = min(os_mean_inh+t2-t1 , length(wtModel));
%     exhal_wtModel = wtModel((os_mean_inh+1):t_end);
%     
%     modelLength = length(inhal_wtModel)+length(exhal_wtModel);
%     warpTimeModel(t_ind+(1:modelLength)) = [inhal_wtModel exhal_wtModel'];

    %%%%%% volume model %%%%%%
    ck = correction(vec(i))-mean_corr;
    if ck < 0
        t_vec = abs(ck) + (1:t2);
        if t_vec(end)>length(volModel)
            volModel = [volModel; zeros(t_vec(end)-length(volModel), 1)];
        end
        volumeModel(t_ind+1:(t_ind+t2)) = volModel(t_vec);
    elseif ck > 0
        t_vec = 1:t2;
        volumeModel(t_ind+1:(t_ind+t2)) = [zeros(abs(ck), 1); volModel(t_vec(1:(t2-ck)))];
    else
        t_vec = 1:t2;
        volumeModel(t_ind+1:(t_ind+t2)) = volModel(t_vec);
    end
            
    %%%%%% time-dilation model %%%%%%
%     tlen = length(tdModel);
%     interp_tdModel = interp1(1:tlen, tdModel, os_mean_inh/t1:os_mean_inh/t1:tlen , 'spline');
%     inhal_tdModel = interp_tdModel(1:os_mean_inh);
%     if  t2 < length(interp_tdModel)
%         exhal_tdModel = interp_tdModel((os_mean_inh+1):t2);
%     else
%         exhal_tdModel = interp_tdModel(os_mean_inh+1:end);
%     end
%     ltdModel = length(inhal_tdModel) + length(exhal_tdModel);
%     timeDilationModel(t_ind+(1:ltdModel)) = [inhal_tdModel exhal_tdModel];


    t_ind =  t_ind + t2;
end

h = hist(find(realTrace),bin/2:bin:(length(realTrace)-bin/2));
% sigma = 15;
% tk = linspace(-sigma*2,sigma*2,(4*sigma)/bin);
% K = exp(-tk.^2/(2*sigma^2));
% K = K/sum(K);
% h = conv(h, K, 'same');
realTrace = reshape(repmat(h, bin, 1), [], 1);

if 0

    gs = buildFrame(55);
    subplot(gs(7))
    
%     plot(sniffTrace/2e4+3,'k');
    plot(realTrace+1,'k');  hold on;
    plot(1:length(timeModel), timeModel, 1:length(warpModel), warpModel-1,...
         1:length(warpTimeModel), warpTimeModel-2, 1:length(volumeModel), volumeModel-3,...
         1:length(timeDilationModel), timeDilationModel-4);
    hold on;
    legend('real data', 'time', 'warp-warp', 'warp-time', 'volume', 'dilation')
    
    t2 = 0;
    for ir = sn_ind
        t1 = ceil(un.sniff_dur(ir));
        t2 = t2+min(t1-rem(t2, bin), 900);
        plot([t2 t2], [-3 2], 'k--')
    end
    hold off
    set(gs(7),  'YLim',        [-4 4],...
                'YTick',       [],...
                'YTickLabel',  [],...
                'box',         'off')
    
    % plot of 4 responses
    subplot(gs(6))
    p = plot((1:length(rt_psth))*bin, rt_psth, ...
             (1:length(rw_psth))*bin, rw_psth, ...
             (1:length(rwt_psth))*bin, rwt_psth, ...
             (1:length(rvol_psth))*bin, rvol_psth, ...
             (1:length(rdil_psth))*bin, rdil_psth);
    set(gca,  'XLim',        [0 500],...
                    'XTick',       [],...
                    'XTickLabel',  [],...
                    'YLim',        [0 max(rt_psth)],...
                    'YTick',       [0 max(rt_psth)],...
                    'YTickLabel',  [0 max(rt_psth)],...
                    'box',         'off')  
                
    [y, x]   = find(rt'==1);
    [y1, x1] = find(rw'==1);
    [y2, x2] = find(rwt'==1);
    [y3, x3] = find(rvol'==1);
    [y4, x4] = find(rdil'==1);


    subplot(gs(1))
    plot(x, y, 'Marker', '.', 'LineStyle', 'none', 'Color', p(1).Color); hold on
    for i = 1:length(inh_ind)
        plot(un.inh_dur(inh_ind(i)), i, 'k.')
%         plot(un.sniff_dur(rind(i)), i, 'k.')
    end
    hold off
    xlim([0, max(x)]);
    legend('time')

    subplot(gs(2))
    plot(x1, y1, 'Marker', '.', 'LineStyle', 'none', 'Color', p(2).Color)
    xlim([0, max(x)]);
    legend('warp')

    subplot(gs(3))
    plot(x2, y2, 'Marker', '.', 'LineStyle', 'none', 'Color', p(3).Color)
    xlim([0, max(x)]);
    legend('wt')
    
    subplot(gs(4))
    plot(x3, y3, 'Marker', '.', 'LineStyle', 'none', 'Color', p(4).Color)
    xlim([0, max(x)]);
    legend('vol')
    
    subplot(gs(5))
    plot(x4, y4, 'Marker', '.', 'LineStyle', 'none', 'Color', p(5).Color)
    xlim([0, max(x)]);
    legend('td')
    

    for i = 1:5
        set(gs(i),  'XLim',        [0 500],...
                    'XTick',       [],...
                    'XTickLabel',  [],...
                    'YLim',        [0 max(y1)+1],...
                    'YTick',       [0 10*floor(max(y1)/10)],...
                    'YTickLabel',  [0 10*floor(max(y1)/10)],...
                    'box',         'off')
    end
end

writeVersionInfo();

% r^2 calculation
var_realTrace = sum(realTrace.^2);
r2_rt   = 1 - sum((realTrace' - timeModel(1:length(realTrace))').^2)/var_realTrace;
r2_rw   = 1 - sum((realTrace' - warpModel(1:length(realTrace))').^2)/var_realTrace;
r2_rwt  = 1 - sum((realTrace' - warpTimeModel(1:length(realTrace))') .^2)/var_realTrace;
r2_vol  = 1 - sum((realTrace' - volumeModel (1:length(realTrace))') .^2)/var_realTrace;
r2_rdil = 1 - sum((realTrace' - timeDilationModel(1:length(realTrace))').^2)/var_realTrace;


t_vec = 1:2e3;
% figure;
% plot(t_vec, [realTrace(t_vec), timeModel(t_vec), volumeModel(t_vec)]);
ll_rt   = logLikelihood(realTrace, timeModel, t_vec);
ll_rw   = logLikelihood(realTrace, warpModel, t_vec);
ll_rwt  = logLikelihood(realTrace, warpTimeModel, t_vec);
ll_vol  = logLikelihood(realTrace, volumeModel, t_vec);
ll_rdil = logLikelihood(realTrace, timeDilationModel, t_vec);


fprintf('r^2 calc: rt: %5.3f,  rw: %5.3f, rwt: %5.3f, rvol: %5.3f, rdil: %5.3f\n',...
                    r2_rt, r2_rw, r2_rwt, r2_vol, r2_rdil);
fprintf('LL  calc: rt: %5.3f,  rw: %5.3f, rwt: %5.3f, rvol: %5.3f, rdil: %5.3f\n',...
                    ll_rt, ll_rw, ll_rwt, ll_vol, ll_rdil);
end

function ll = logLikelihood(realData, model, t_vec)

    model = model*.99+.005;
    ll    = sum(realData(t_vec).*log(model(t_vec)) - model(t_vec));
end

function [Lvec, r2vec] = bestL_Estimation(un, md, sn, bin, maxPerc)
%
%
%
my = my_tools();

% looking for the first sniff cycles for a given odor and given concentration
rind = find (un.valve == md.valve & un.conc == md.conc & un.odorSniffNum == 1);

[inh_dur, ind]  = sort(un.inh_dur(rind));
inh_ind = rind(ind);

Lvec = zeros(11, 1);
r2vec = zeros(11, 1);
perc = 2;

for ip = (0:perc:maxPerc)/perc
    % rectification of negative pressure part
    % integral calculation for each sniff
    % search for time shift
    
    correction = timeShiftByVolume(sn, inh_ind, inh_dur, perc*ip);
    mean_corr = round(mean(correction));
    
    % spike times realignment according to different models
    [rt, ~, ~, rvol, ~]  = spikeRealignment(un, rind, correction);

    %%%%%% time model %%%%%%
    [rt_psth, ~]     = my.nan_spike_count(rt', bin); % in previous version I used nanpsth
    tModel = reshape(repmat(rt_psth, bin, 1),[], 1);
    
    %%%%%% volume model %%%%%%
    [rvol_psth, ~]    = my.nan_spike_count(rvol', bin);
    volModel = reshape(repmat(rvol_psth, bin, 1), [], 1);
    
    
    realTraceDuration   = sum(un.sniff_dur(rind));
    sniffTrace          = zeros(realTraceDuration, 1);
    realTrace           = zeros(realTraceDuration, 1);
    timeModel           = zeros(realTraceDuration, 1);
    volumeModel         = zeros(realTraceDuration, 1);
    
    
    vec = 1:length(inh_ind);
    sn_ind = inh_ind;
    t_ind = 0;
    i = 0;
    for ir = sn_ind
        i = i+1;
        t2 = ceil(un.sniff_dur(ir));
        t2 = min(t2-rem(t2, bin), 900);
        
        %%%%%% sniff trace %%%%%%
        sniffTrace(t_ind+(1:300)) = sn.waveform(ir, :);
        
        %%%%%% real trace %%%%%%
        realTrace(t_ind+1:(t_ind+t2)) = full(un.rt(1:t2, ir))';
        
        %%%%%% time model %%%%%%
        timeModel(t_ind+1:(t_ind+t2)) = tModel(1:t2);
        
        %%%%%% volume model %%%%%%
        try
            ck = correction(vec(i))-mean_corr;
            if ck < 0
                t_vec = abs(ck) + (1:t2);
                if t_vec(end)>length(volModel)
                    volModel = [volModel; zeros(t_vec(end)-length(volModel), 1)];
                end
                volumeModel(t_ind+1:(t_ind+t2)) = volModel(t_vec);
            elseif ck > 0
                t_vec = 1:t2;
                volumeModel(t_ind+1:(t_ind+t2)) = [zeros(abs(ck), 1); volModel(t_vec(1:(t2-ck)))];
            else
                t_vec = 1:t2;
                volumeModel(t_ind+1:(t_ind+t2)) = volModel(t_vec);
            end
            
        catch
            fprintf('Problem in volume model, trial # %03d\n', i)
            keyboard
        end
        
        t_ind =  t_ind + t2;
    end
    
    h = hist(find(realTrace),bin/2:bin:(length(realTrace)-bin/2));
    realTrace = reshape(repmat(h, bin, 1), [], 1);
    
    if 0
        
        gs = buildFrame(55);
        subplot(gs(7))
        
        plot(sniffTrace/2e4+3,'k'); hold on;
        plot(realTrace/2+1,'k');
        plot(1:length(timeModel), timeModel, 1:length(volumeModel), volumeModel-1);
        
        hold on;
        legend('real data', 'time', 'volume')
        
        t2 = 0;
        for ir = sn_ind
            t1 = ceil(un.sniff_dur(ir));
            t2 = t2+min(t1-rem(t2, bin), 900);
            plot([t2 t2], [-3 2], 'k--')
        end
        hold off
        set(gs(7),  'YLim',        [-4 4],...
            'YTick',       [],...
            'YTickLabel',  [],...
            'box',         'off')
        
        % plot of 4 responses
        subplot(gs(6))
        p = plot((1:length(rt_psth))*bin, rt_psth, ...
            (1:length(rvol_psth))*bin, rvol_psth);
        set(gca,  'XLim',        [0 500],...
            'XTick',       [],...
            'XTickLabel',  [],...
            'YLim',        [0 max(rt_psth)],...
            'YTick',       [0 max(rt_psth)],...
            'YTickLabel',  [0 max(rt_psth)],...
            'box',         'off')
        
        [y, x]   = find(rt'==1);
        [y3, x3] = find(rvol'==1);
        
        subplot(gs(1))
        plot(x, y, 'Marker', '.', 'LineStyle', 'none', 'Color', p(1).Color); hold on
        for i = 1:length(inh_ind)
            plot(un.inh_dur(inh_ind(i)), i, 'k.')
            %         plot(un.sniff_dur(rind(i)), i, 'k.')
        end
        hold off
        xlim([0, max(x)]);
        legend('time')
        
        subplot(gs(4))
        plot(x3, y3, 'Marker', '.', 'LineStyle', 'none', 'Color', p(4).Color)
        xlim([0, max(x)]);
        legend('vol')
        
        for i = 1:5
            set(gs(i),  'XLim',        [0 500],...
                'XTick',       [],...
                'XTickLabel',  [],...
                'YLim',        [0 max(y1)+1],...
                'YTick',       [0 10*floor(max(y1)/10)],...
                'YTickLabel',  [0 10*floor(max(y1)/10)],...
                'box',         'off')
        end
    end
    
    % r^2 calculation
    var_realTrace = sum(realTrace.^2);
    r2_rt   = 1 - sum((realTrace' - timeModel(1:length(realTrace))').^2)/var_realTrace;
    r2_vol  = 1 - sum((realTrace' - volumeModel (1:length(realTrace))') .^2)/var_realTrace;
    r2vec(ip+1) = r2_vol;
    Lvec (ip+1) = ip*perc;
%     fprintf('r^2 calc: rt: %5.3f,  rvol: %5.3f\n', r2_rt,  r2_vol);
end
end

function [rt, rw, rwt, rvol, rdil]  = spikeRealignment(un, rind, correction)

[~, ind]  = sort(un.inh_dur(rind));
inh_ind = rind(ind);

%calculation of mean parameters
mean_inh  = round(mean(un.inh_dur));
mean_sn   = round(mean(un.sniff_dur));
mean_exh  = mean_sn - mean_inh;

% os - is for odor specific
os_mean_inh  = round(mean(un.inh_dur(rind)));
os_mean_sn  = round(mean(un.sniff_dur(rind)));
os_mean_exh = os_mean_sn - os_mean_inh;

rt   = zeros(max(ceil(un.sniff_dur(rind))), length(rind));
rw   = zeros(mean_inh + mean_exh, length(rind));
rwt  = zeros(mean_inh + mean_exh, length(rind));
rvol = zeros(max(ceil(un.sniff_dur(rind))), length(rind));
rdil = zeros(mean_inh + mean_exh, length(rind));

mean_corr = round(mean(correction));

k = 0;  
yind = 0;      
    % spike times realignment according to different models
    for ir = inh_ind
        k = k+1;
        yind = yind+1;

        t1 = round(un.inh_dur(ir));
        t2 = ceil(un.sniff_dur(ir));

        if t2<=size(un.rt, 1)

%             inh_spikes = find(full(un.rt(1:t1, ir))==1);
%             exh_spikes = find(full(un.rt((t1+1):t2, ir))==1);
            all_sp = find(full(un.rt(1:t2, ir))==1);

            % 'time'
            rt(1:t2, yind) = un.rt(1:t2, ir);
            rt((t2+1):end, yind) = NaN;

%             % 'warp-warp'
%             rw(round(inh_spikes/t1*os_mean_inh), yind) = 1;
%             rw(round(exh_spikes/(t2-t1)*os_mean_exh + os_mean_inh), yind) = 1;

%             % 'warp-time'
%             rwt(round(    inh_spikes/t1*(os_mean_inh)), yind) = 1;
%             rwt(round(exh_spikes + round(os_mean_inh)), yind) = 1;
%             rwt((os_mean_inh+(t2-t1)+1):end, yind) = NaN;

            % 'volume' : aligment of the spike train by inhalation volume
            sp = round(all_sp) - correction(k) + mean_corr;
            rvol(sp(sp>0), yind) = 1;
            rvol((t2+1):end, yind) = NaN;
            % rvol((os_mean_inh+(t2-t1)+1):end, yind) = NaN;

%             % 'time dilation'
%             rdil_spikes = round(all_sp/t1*mean_inh);
%             rdil(rdil_spikes, yind) = 1;
%             rdil((max(rdil_spikes)+1):end, yind) = NaN;
        end
    end

end

function correction = timeShiftByVolume(sn, inh_ind, inh_dur, perc)

% rectification of negative pressure part
sn.waveform(sn.waveform<0)=0;
% integral calculation for each sniff
maxVolume = zeros(1, length(inh_dur));
for k = 1:length(inh_dur)
    try
        volume = cumsum(int32(sn.waveform(inh_ind(k), :)));
    catch
        keyboard
    end
    maxVolume(k) = max(volume);
end


% search for time shift
                                % mv = sort(maxVolume);
tshift = mean(maxVolume)/20;    % mv(end-3)/20;

correction = zeros(1, length(inh_dur));
for k = 1:length(inh_dur)
    try
        volume = cumsum(int32(sn.waveform(inh_ind(k), :)));
        if max(volume) < perc/5*tshift
            correction(k) = round(inh_dur(k));
        else
            correction(k) = find(volume>perc/5*tshift, 1);
        end
    catch
        fprintf('Problem with volume calculation\n')
        keyboard
    end
end
end

function gs = buildFrame(fig)
    f = figure(fig);
    clf;
    set(f, 'Color', 'w')

    bvec = [0.8:-0.14:0.24, 0.07];
    hvec = [0.13*ones(1, 6) 0.2] ;
    width = [0.25 0.5, 0.12];
    lvec = [0.05, 0.32, 0.85];
    gs = zeros(1, 16);
    for i = 1:6
        gs(i)  = subplot('position', [lvec(1) bvec(i) width(1) hvec(i)]);
    end
    gs(7)  = subplot('position', [lvec(2) 0.07 0.6 0.7]);
   
end

function writeVersionInfo()
% writing info about given figure creation
fullPath = [mfilename('fullpath') '.m'];
partPath = strsplit(fullPath, 'box');
DirInfo = dir(fullPath);
date = DirInfo.date;
str = sprintf('Was generated using %s from %s', char(partPath(2)), date);
set(gcf, 'NumberTitle', 'off')
set(gcf, 'Name', str)
end

function [m_psth, std_psth] = meanAndStd(sp, bin, range)
    
    range(2) = min(range(2), size(sp, 2));
    [m_psth, std_psth] = psth_and_error(sp(:, range(1):range(2))', bin);

    function [mean_y, std_y] = psth_and_error(rast, binsize)
        
        np  = size(rast, 2);
        
        if rem(np,binsize)>0
            np = np - rem(np,binsize);
        end
        
        [y, x] = find(rast);
        
        for j = 1:max(y)
            psth(j, :) = histcounts(x(y==j), 0:binsize:np)/binsize*1e3;
        end
        mean_y = mean(psth, 1);
        std_y  = std(psth, [], 1);
        
    end
end

function h = plotRasterAndPSTH(resp, bin, lat)
    my = my_tools();
    
    rvec = (lat-30):(lat+20);
    
%     [x, y] = my.raster2plot(resp);
%     mean_psth = my.psth(resp, bin);

    hist2D = my.rast2hist(resp(:, rvec), bin);
    h = mean(entropy(hist2D));

%     plot(x, y, '.', 'Color', [0, 0.7, 0]); hold on;
%     plot(bin*(1:length(mean_psth)), mean_psth, 'k'); 
% %     plot(bin*(1:length(ci_low)), [ci_low ci_high], 'k--'); 
% %     plot(bin*(1:length(ci_low)), 50*h, 'k--'); 
%     hold off;
% 
%     legend(sprintf('%4.2f',mean(h)))
%     set(gca, 'box', 'off', 'YTick', []);
%     ylim([0 max(y)+1]); xlim([0 600]);
%     drawnow
end

function H = entropy(x)

m = size(x, 2); %number of bins
H = zeros(1,m);

for col = 1:m
    % list of unique values in each bin
    listOfVal = unique(x(:,col));
    freq = zeros(size(listOfVal));
	
    % Calculate sample frequencies
    for num = 1:length(listOfVal)
        freq(num) = sum(x(:,col) == listOfVal(num));
    end
	
    % Calculate sample class probabilities
    P = freq / sum(freq);
    H(col) = -sum(P .* log2(P));
end

end

function diagHist(x,y, color, shift)
    xvec = -30:2.5:60;
    xshift = max(xvec) + min(xvec);

    h1 = fliplr(hist((y-x)/sqrt(2), xvec/sqrt(2)));
    [xx, yy] = stairs((xvec - xshift)/sqrt(2), h1);
    P = [cos(-pi/4) -sin(-pi/4); sin(-pi/4)   cos(-pi/4)];
    R = P*[xx'; yy'];

    plot( R(1,:)+shift, R(2,:)+shift, 'Color', color)

end

function offDiag = diagBoxPlot(HL, H, L, col, ind, shift, respType)
    if ~isempty(L)
        y1 = abs(HL-H);
        x1 = abs(H-L);
    else
        y1 = HL';
        x1 = H';
    end
%     sort((y1-x1)/sqrt(2))
%     mean(sort((y1-x1)/sqrt(2)))

%     meann = mean(xx1);
%     stdd = std(xx1);
%     I = bsxfun(@gt, abs(bsxfun(@minus, xx1, meann)), 2*stdd);
%     out = find(I);

    xx1 = -prctile((y1-x1)/sqrt(2),[25 50 75],2);
    yy1 = [ind ind ind];
    
    P = [cos(-pi/4)  -sin(-pi/4); sin(-pi/4)   cos(-pi/4)];
    R = P*[xx1; yy1];

    plot( R(1,:)+shift, R(2,:)+shift, 'Color', col, 'LineWidth', 2)
    plot( R(1,2)+shift, R(2,2)+shift, 'Color', col, 'Marker', '.', 'MarkerSize', 20)
    
    offDiag = (y1-x1)/sqrt(2);
    [~, p_val] = ttest(offDiag);
    fprintf('ttest %s sp. count, p-value: %0.3g\n', respType, p_val);
    [ p_val, ~] = signrank(offDiag);
    fprintf('Wilcoxon signed rank test %s sp. count, p-value: %0.3g\n', respType, p_val);
end
