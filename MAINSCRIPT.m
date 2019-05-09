clear,close all;

% Loop parameters
numruns = 100;
numtimeslots = 100;
nlist = [1 9 19 29 39 49];
simavg = zeros(length(nlist),1);
simavg_d = zeros(length(nlist),1);

% System parameters
SINR_thresh = 5;
Rd = 1; % desired communication rate (bits/channel use)
% N = 1;          % number of chips per symbol
N_list = [1 4 8 16];
spc = 4;        % samples per chip
Ns = 20;        % number of symbols per packet
pulseshape = 'rect';    % pulse shape (either 'rect' or 'rrc')
x = 10;         % 'x' parameter, see G. Kuperman paper on MB-URAM
M = 4;          % number of antennas
L = 2;        % number of beams each node can form
m = [0:M-1].';  % antenna indices
Fc = 300e6;     % carrier frequency [Hz]
c = physconst('LightSpeed');    % speed of light [m/sec]
lambda = c/Fc;  % wavelength [m]
d = lambda/2;   % antenna spacing [m]
Tc = 1/(4e3);   % chip period [sec]
Ts = Tc/spc;    % sampling rate [Hz]
Pnoise = physconst('Boltzman')*300/Tc; % noise power [watts]
SNR_dB = 5;     % minimum SNR [dB]
pl = 2;  % power loss coefficient: pl = 2 for free space, pl = 3 for underwater
ARRIVAL_EVENT = 1;
TERMINATION_EVENT = 0;
radius = 100;  % radius of simulation circle [meters]
allowed_attempts = 100;
beam_opt = 'independent';   % either 'independent' or 'identical'
s = string(datetime('now','Format','dMMMy_HHmm'));

% Build pulse shaping filter
switch pulseshape
    case 'rrc'
        p = rcosdesign(rolloff,span,spc,shape);
    case 'rect'
        p = rectpulse(1,spc);
end

% Loop over varying chips/symbol
for chips = 1:length(N_list)
    N = N_list(chips);
    fprintf('N = %d: ', N);
    spp = spc*N*Ns;
    save_file = sprintf('%s_r%d_L%d_M%d_N%d', s, pl, L, M, N);

    % Loop over number of nodes list
    for nni = 1:length(nlist)
        fprintf('%d ',nlist(nni));
        nn = nlist(nni);

        % Loop for number of runs
        for run = 1:numruns    
            Rd = 1;
            dFlag = 0;
            txAttempts = 0;
            txSuccess = 0;
            txSuccess_d = 0;
            MBURAM_success = 0;
            MBURAM_att = 0;
            SINR_MF_hist = [];
            SINR_D_hist = [];
            Rmf_hist = [];
            Rd_hist = [];

            MF_attempts = 0;
            DEC_attempts = 0;

            numerr_d = 0;
            numerr_mf = 0;

            % Initialization of parameters and system state
            [ positions, angles, ulaAng ] = createNodeLayout(nn, radius);
            n = size(positions,1)-1;
            Ptx = getTxPower(positions,SNR_dB,Pnoise,pl);

            nb = min(L, n);
            codeInd = repelem(1:(n+1), nb);
            code = makeCodes(N, nb, n+1);

            packetsRxd = 0;
            packetsRxd_D = 0;
            activeLinks = zeros(n+1,1);
            rxingLink = zeros(n+1,1);
            activeChannels = zeros(n+1,1);
            rxingChannels = zeros(n+1,1);
            t = 0;
            idletimes = (randi([0, spp],n+1,1));
            for ii = 1:length(idletimes)
               r = norm(positions(1,:)-positions(ii,:));
               tau = r/physconst('LightSpeed');
               del = round(tau/Ts);
               idletimes(ii) = idletimes(ii) + del; 
            end

            state = [];
            state(idletimes == 0) = 1;
            state(idletimes ~= 0) = 0;
            typeList = [];
            txnodeList = [];
            tstartList = [];
            taulist = cell(n+1,1);
            txs = find(state);
            for txnode = txs
                tstart = idletimes(txnode);
                tend = tstart + spp;
                typeList = [ typeList ARRIVAL_EVENT ];
                txnodeList = [ txnodeList txnode ];
                tstartList = [ tstartList tstart ];
                typeList = [ typeList TERMINATION_EVENT ];
                txnodeList = [ txnodeList txnode ];
                tstartList = [ tstartList tend ];
            end

            rxs = find(~state);
            for rxnode = rxs
                tstart = idletimes(rxnode);
                tend = tstart + spp;
                typeList = [ typeList ARRIVAL_EVENT ];
                txnodeList = [ txnodeList rxnode ];
                tstartList = [ tstartList tstart ];
                typeList = [ typeList TERMINATION_EVENT ];
                txnodeList = [ txnodeList rxnode ];
                tstartList = [ tstartList tend ];
            end

            [ tstartList, ind ] = sort( tstartList );
            typeList = typeList( ind );
            txnodeList = txnodeList( ind );

            % Main simulation loop
            while txAttempts < allowed_attempts
                eventType = typeList(1);
                txnode = txnodeList(1);
                t = tstartList(1);
                typeList = typeList(2:end);
                txnodeList = txnodeList(2:end);
                tstartList = tstartList(2:end);

                state(idletimes > t) = 0;
                state(idletimes <= t) = 1;

                % if node 1 is due to transmit, then all of the links it was receiving
                % on are losts
                if txnode == 1
                    activeLinks = zeros(n+1,1);
                    activeChannels = zeros(n+1,1);
                end

                % only considering arrivals at node 1 for simplicity
                if isequal(eventType,ARRIVAL_EVENT)
                    rxingLink(txnode) = 1;
                    picks = chooseLNodes(L, [1 n+1], txnode);
                    MBURAM_att = MBURAM_att + 1;

                    if state(1) ~= 1 && ~isempty(picks == 1)
                        taulist{txnode} = idletimes - t;
                        activeLinks(txnode) = 1;
                    end
                elseif isequal(eventType,TERMINATION_EVENT)
                    if state(1) ~= 1
                        MBURAM_success = MBURAM_success + 1;
                        if activeLinks(txnode)

                            txAttempts = txAttempts + 1;

                            tau = taulist{txnode};
                            users = find(tau+1 <= spp);
                            uusers = users;
                            ind = find(tau(users) < 0);
                            for i = 1:length(ind)
                                if idletimes(users(ind(i)))-t+spp <= spc*N*Ns
                                    tau(end+1) = idletimes(users(ind(i)))-t+spp;
                                    uusers(end+1) = uusers(ind(i));
                                end
                            end

                            r = sqrt((positions(1,1)-positions(:,1)).^2+(positions(1,2)-positions(:,2)).^2);
                            Prx_all = Ptx./r.^pl; Prx = Prx_all(uusers);                
                            g = zeros(length(uusers),1);

                            ula = ulaAng; ang = angles;
                            ind = find(ang > 180);
                            ang(ind) = ang(ind) - 360;
                            phi_rx0 = getaoa(txnode,1,ula,ang);
                            tau_rx = m*d*sin(deg2rad(phi_rx0))/c;
                            cmf_rx = exp(-1j*2*pi*Fc*tau_rx)/sqrt(M);

                            % g is a vector with the receive beam gain for each
                            % of the transmitting users
                            for i = 1:length(uusers)
                                u = uusers(i);
                                if abs(ang(1, u) - ang(1, txnode)) <= 90
                                    phi_rx = getaoa(u,1,ula,angles);
                                    tau_rx = m*d*sin(deg2rad(phi_rx))/c;
                                    v_rx = exp(-1j*2*pi*Fc*tau_rx)/sqrt(M);
                                    F_rx = cmf_rx'*v_rx;
                                    g(uusers==u) = abs(F_rx)^2;
                                else
                                    g(uusers==u) = 0;
                                end
                            end

                            picks = [];
                            gg = [];
                            % For each transmitting node
                            for uu = 1:length(uusers)
                                u = uusers(uu);

                                % Choose L nodes to transmit to.
                                picks(:,uu) = chooseLNodes(L,[1 n+1],u)';
                                if n+1 > L
                                    if activeLinks(u) == 1
                                        picks(1,uu) = 1;                               
                                    else
                                        while picks(1,uu) == 1
                                            picks(:,uu) = chooseLNodes(L,[1 n+1],u)';
                                        end
                                    end
                                else
                                    picks(:,uu) = chooseLNodes(L,[1 n+1],u)';
                                end

                                phi_tx = getaoa(1,u,ula,angles);
                                tau_tx = m*d*sin(deg2rad(phi_tx))/c;
                                v_tx = exp(-1j*2*pi*Fc*tau_tx)/sqrt(M);
                                cmf_tx = v_tx;
                                F_tx = cmf_tx'*v_tx;

                                % gg(rx, tx) corresponds to transmit beam gain
                                % in direction of central node. For L = 1, all
                                % gains will be = 1 where TX steers directly at
                                % central node. For L > 1, there will be more
                                % rows in gg corresponding to "out-of-beam"
                                % interference
                                for uuu = 1:size(picks,1)
                                    jj = picks(uuu,uu);
                                    if abs(ang(u, jj) - ang(u, 1)) <= 90
                                        phi_txjj = getaoa(jj,u,ula,angles);
                                        tau_txjj = m*d*sin(deg2rad(phi_txjj))/c;
                                        v_txjj = exp(-1j*2*pi*Fc*tau_txjj)/sqrt(M);
                                        F_txjj = cmf_tx'*v_txjj;
                                        gg(uuu,uu) = abs(F_txjj)^2; 
                                    else
                                        gg(uuu,uu) =0 ;
                                    end
                                end
                            end
                            gtmp = repmat(g', [size(gg,1) 1]);
                            g = (gg.*gtmp);
                            g = g(:);   % at this point, g is a vector of length(users)
                            Prxx_uu = [];
                            for ii = 1:length(uusers) % now make g a vector of length(uusers)
                                Prxx_uu(ii) = Prx(uusers(ii) == users);
                            end
                            Prxx = repelem(Prxx_uu(:),size(picks,1)); % account for possible multiple beams
                            Prxx = Prxx(:).*g;   % apply beam forming power gain

                            % full A matrix with everything in it
                            A = diag(repmat(sqrt(Prxx(:)),[Ns 1]));

                            % K is the total number of signals = number of
                            % active TXs times L (the number of beams each is sending).
                            % One of those signals is the SOI, the rest are
                            % either signals intended for a different RX or
                            % signals coming from other TXs
                            K = length(g); 

                            totalbeams = users;
                            txind = (find(totalbeams == txnode)-1)*size(picks,1)+1;
                            codeb = code(:, ismember(codeInd, users));
                            taubeams = repelem(tau(users), size(picks, 1));
                            idletimes_b = repelem(idletimes(users), size(picks, 1));
                            [S, ~] = getSMatrix(spc, Ns, taubeams, idletimes_b-t+spp, codeb, p, repelem(users,size(picks,1)));
                            S = S(abs(min(taubeams))+1:abs(min(taubeams))+spp,:); % trim out only SOI's packet duration

                            % delete any column with all zeros, corresponding
                            % to a transmitted symbol that doesn't overlap at
                            % all in time, happens because of partial overlap
                            % between packets
                            cols2del = find(sum(abs(S)) == 0);
                            S(:, cols2del) = [];

                            % This piece of code goes through and provides labels to easily read the user and beam 
                            % that each column belongs to. Then soi_ind and
                            % int_ind are the indices of the SOI and all the
                            % interference signals present
                            uusers_new = repmat(repelem(uusers,nb,1),[Ns,1]);
                            uusers_new(cols2del) = [];
                            uusers_new_cell = [];
                            for ii = 1:nb:length(uusers_new)
                                for jj = 0:nb-1
                                    uusers_new_cell{ii+jj} = sprintf('U%d_B%d',uusers_new(ii+jj),jj+1);
                                end
                            end
                            soi_ind = [];
                            int_ind = [];
                            for ii = 1:length(uusers_new_cell)
                                if isequal(uusers_new_cell{ii}, sprintf('U%d_B1', txnode))
                                    soi_ind(end+1) = ii;
                                else
                                    int_ind(end+1) = ii;
                                end
                            end

                            % Create the A matrix for the SOI and interference
                            Ad = diag(A);
                            Ad(cols2del) = [];
                            Asoi = diag(Ad(soi_ind));
                            Aint = diag(Ad(int_ind));

                            % Create the SOI and interference BPSK symbols
                            b = 2*randi([0, 1],size(S, 2),1)-1;
                            bsoi = b(soi_ind);
                            bint = b(int_ind);

                            % Split the S matrix into SOI and interference
                            % pieces
                            Ssoi = S(:, soi_ind);
                            Sint = S(:, int_ind);

                            % Total received DT signal is: ysoi + yint + w
                            ysoi = Ssoi*Asoi*bsoi;
                            yint = Sint*Aint*bint;
                            w = sqrt(Pnoise)*randn(size(ysoi));

                            % Matched filter outputs
                            rsoi = Ssoi'*ysoi;
                            rint = Ssoi'*yint;
                            rw = Ssoi'*w;
                            % trick to be able to run less realizations since noise power at MF input = noise power at MF output
                            rw = repmat(sqrt(Pnoise), size(rw)); 

                            % if there is no interference just make sure MF
                            % outputs zeros for interference component
                            if isempty(rint)
                                rint = zeros(size(rsoi));
                            end

                            % Record the SINR and Rate information
                            SINR_MF = (rsoi.^2) ./ (rint.^2 + rw.^2);        
                            SINR_MF_hist(:, end+1) = SINR_MF;                     
                            Rmf = .5*log2(1+SINR_MF);
                            Rmf_hist(:, end+1) = Rmf;

                            % Record the number of symbol errors (= bit errors
                            % for BPSK)
                            numerr_mf = numerr_mf + length(find(Rmf < 1));
                            MF_attempts = MF_attempts + 1;

                            % If Rate is greater than desired rate --> success
                            if Rmf >= Rd
                                txSuccess = txSuccess + 1;
                                txSuccess_d = txSuccess_d + 1;
                                packetsRxd = packetsRxd + 1;
                                packetsRxd_D = packetsRxd_D + 1;
                                Rd_hist(:, end+1) = Rmf;
                                SINR_D_hist(:, end+1) = SINR_MF;
                                numerr_d = numerr_d + length(find(Rmf < 1));
                            else    % MF Rx drops the packet, now try MUD

                                DEC_attempts = DEC_attempts + 1;

                                % This step uses checks each node to determine
                                % which nodes are causing "strong enough"
                                % interference. We consider a node's
                                % interference to be "strong enough" if it
                                % causes the SINR at the Rx to fall within 3 dB
                                % of the minimum required for a packet drop,
                                % 0.5*log2(1+SINR) = 1 bit/channel use --> SINR
                                % = 3 = 4.77 dB  
                                Psoi = Asoi(1, 1)^2;
                                mud_ind = find(Psoi./(diag(Aint).^2 + Pnoise) < SINR_thresh);

                                % Gather columns to include in the MUD process
                                % SOI is always included in the MUD
                                Amud_int = diag(Aint);
                                Amud_int = diag(Amud_int(mud_ind));
                                Smud_int = Sint(:, mud_ind);
                                bmud_int = bint(mud_ind);
                                Smud = [Ssoi, Smud_int];

                                % Set B = the left over columns which are not
                                % included in the MUD process and treated as
                                % noise
                                A_B = diag(Aint);
                                A_B(mud_ind) = [];
                                A_B = diag(A_B);
                                S_B = Sint;
                                S_B(:, mud_ind) = [];
                                b_B = bint;
                                b_B(mud_ind) = [];
                                % in case all columns are included, just make
                                % sure that set B contributes nothing otherwise
                                % generate the signal
                                if isempty(b_B)
                                    y_B = zeros(size(Smud,1),1);
                                else
                                    y_B = S_B*A_B*b_B;
                                end

                                % correlation matrix (MF outputs)
                                R = Smud'*Smud;
                                try
                                    F = pinv(R); % decorrelator filter bank
%                                     F_dt = pinv(Smud'*Smud)*Smud'; % this decorrelator operates directly on DT signals (ysoi, ymud_int, etc.)
                                catch
                                    dFlag = 1; % if pinv() never converges (very rarely happens)
                                end
                                if ~dFlag                                  
                                    ysoi = Ssoi*Asoi*bsoi;                                
                                    ymud_int = Smud_int*Amud_int*bmud_int;
                                    if isempty(ymud_int)
                                        ymud_int = zeros(size(ysoi));
                                    end

                                    % Perform matched filtering
                                    rsoi = Smud'*ysoi;
                                    rint = Smud'*ymud_int;
                                    rw = repmat(sqrt(Pnoise), size(rint));
                                    r_B = Smud'*y_B;

                                    % Perform decorrelator and get entries
                                    % corresponding to the signal of interest
                                    xsoi = F*rsoi; xsoi = xsoi(1:Ns);
                                    xint = F*rint; xint = xint(1:Ns); % should be very small values if the MUD was able to do its job, otherwise R was likely ill conditioned
                                    xw = F*rw; xw = xw(1:Ns);                         
                                    x_B = F*r_B; x_B = x_B(1:Ns);
                                end
                                if ~dFlag
                                    SINR_D = xsoi.^2./(x_B.^2 + xw.^2 + xint.^2);
                                    SINR_D_hist(:, end+1) = SINR_D;                                
                                    Rdec = .5*log2(1+SINR_D);
                                    numerr_d = numerr_d + length(find(Rdec < 1));
                                    Rd_hist(:, end+1) = Rdec;
                                    if Rdec >= Rd
                                        txSuccess_d = txSuccess_d + 1;
                                        packetsRxd_D = packetsRxd_D + 1;
                                    end
                                else
                                    Rdec = 0; % just to ensure the packet doesn't get counted
                                end
                            end
                        end
                    end

                    % disconnect the link
                    activeLinks(txnode) = 0;
                    rxingLink(txnode) = 0;

                    % get a new idle time and add event to the queue
                    idletimes(txnode) = t + randsrc(1, 1, [0, spp*x; (x-1)/x, 1/x]);
                    tstart = idletimes(txnode);
                    tend = tstart + spp;
                    typeList = [typeList, ARRIVAL_EVENT];
                    txnodeList = [txnodeList, txnode];
                    tstartList = [tstartList, tstart];
                    typeList = [typeList, TERMINATION_EVENT];
                    txnodeList = [txnodeList, txnode];
                    tstartList = [tstartList, tend];
                    [tstartList, ind] = sort(tstartList);
                    txnodeList = txnodeList(ind);
                    typeList = typeList(ind);
                end
            end

            % Generate data points for the run
            sim(run) = spc*Ns*packetsRxd/(Tc*t);
            sim_d(run) = spc*Ns*packetsRxd_D/(Tc*t);
            sim_attempts(run) = txAttempts;
            sim_success(run) = txSuccess;
            sim_success_d(run) = txSuccess_d;

            sim_mburam(run) = MBURAM_success;
            sim_mburamatt(run) = MBURAM_att;

            sim_numerr_d(run) = numerr_d/(Ns*DEC_attempts);
            sim_numerr_mf(run) = numerr_mf/(Ns*MF_attempts);

            sim_sinr_mf_hist{run} = SINR_MF_hist(:);
            sim_sinr_d_hist{run} = SINR_D_hist(:);
            sim_rate_mf_hist{run} = Rmf_hist(:);
            sim_rate_d_hist{run} = Rd_hist(:);
        end

        % Generate data points for the number of neighbors, n
        simavg(nni) = mean(sim);
        simavg_d(nni) = mean(sim_d);
        avg_attempts(nni) = mean(sim_attempts);
        avg_success(nni) = mean(sim_success);
        avg_success_d(nni) = mean(sim_success_d);

        avg_mburam(nni) = mean(sim_mburam);
        avg_mburamsim(nni) = mean(sim_mburamatt);

        avg_err_d(nni) = mean(sim_numerr_d);
        avg_err_mf(nni) = mean(sim_numerr_mf); 

        avg_sinr_mf{nni, :} = sim_sinr_mf_hist;
        avg_sinr_d{nni, :} = sim_sinr_d_hist;
        avg_rate_mf{nni, :} = sim_rate_mf_hist;
        avg_rate_d{nni, :} = sim_rate_d_hist;

    end

    fprintf('\n');
    save(save_file);

end