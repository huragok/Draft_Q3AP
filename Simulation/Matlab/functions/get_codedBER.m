function codedBER = get_codedBER(constellation, map, mu_h, sigma2_h, sigma2_v, max_frame, iter_max, coding_rate, nldpc, seed)
%   codedBER = get_codedBER(constellation, map, mu_h, sigma2_h, sigma2_v, max_frame, iter_max, coding_rate, nldpc, seed)
%   Evaluate the LDPC-coded BER of a specific MoDiv mapping design after
%   the retransmission
% _________________________________________________________________________
%	Inputs:
%       constellation:	Q-by-1 vector, the modulated constellations
%       map:            2-by-Q vector, each row is a permutation of 1 : Q
%                       The first row is the mapping from the source to the
%                       destination and the second reos is from the relay
%                       to the destination.
%       mu_h:           3-by-1 vector, the mean value of the Rician channels
%                       (LOS component)
%       sigma2_h:       3-by-1 vector, the variance of the Rician channels
%                       (fading component)
%       sigma2_v:       scalar, the variance of the received AWGN noise
%                       at the destination
%       max_frame:      Scalar, number of LDPC frames in simulation.
%       iter_max:       Sclar, maximum iteration time within the iterative 
%                       receiver.
%       coding_rate:    coding rate of LDPC, {1/2,2/3,3/4,5/6}.
%       nldpc:          Scalar, bit length after channel coding, 
%                       mod(nldpc,24)=0.
%       seed:           Scalar, seed for the random number generator
%	Outputs:
%		codedBER:		Scalar, the encoded BER by Chase combining the M
%                       transmissions
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 06/01/2015
% Codename: Dunkirk
% _________________________________________________________________________
% The design of this function is based on MIMO_BER_sim function by  
% Weiliang Zeng, 03/04/10.
% _________________________________________________________________________

Q = size(map, 2);
Nbps = round(log2(Q)); % Number of bit per symbol
rng(seed);

max_bit_error = 1000; % Count up to this number of bit error we stop the simulation since BER can be measured accurately enough at this point

% Config LDPC
ind = 0;
max_iterations = 30;
decoder_type = 0;

[H_rows, H_cols, P_matrix] = InitializeWiMaxLDPC(coding_rate, nldpc, ind);
bits_per_frame = length(H_cols) - length(P_matrix);

%generate bit vectors and mapped symbols used in MAP detector
bit_mat = (dec2bin(0 : Q - 1) > '0') + 0;
bit_mat_anti = 1 - 2 * bit_mat; % antipodal matrix, logic 1 is mapped to 1, and logic 0 is mapped to -1
sym_mod_mat = constellation([1 : Q; map]).';
        
% Start the transmission frame by frame
error_all = zeros(1, max_frame);
for i = 1 : max_frame
	if mod(i, 5)==0 % Print a 'x' for every 5 frames
		fprintf('x');
	end
	if mod(i, 400)==0 % Print enter for every 400 frames
		fprintf('\n');
	end

	% source, coding, random interlever and modulation
    data = round(rand(1, bits_per_frame)); % Randomly generated information-bearing bit, 1-by-bits_per_frame
    codewordTemp = LdpcEncode(data, H_rows, P_matrix ); % LDPC encoded codeword, 1-by-nldpc
    codeword_invTemp = randintrlv(codewordTemp, 0);  % Interleaved codeword, 1-by-nldpc
    coded_index = bit2idx(codeword_invTemp, Nbps);

    transmit_Mod = zeros(3, ceil(nldpc / Nbps)); % Modulated symbol for each (re)transmission
    transmit_Mod(1, :) = constellation(coded_index); % Modulation and normalization
    transmit_Mod(2, :) = constellation(map(1, coded_index));
    transmit_Mod(3, :) = constellation(map(2, coded_index));
	
	% Lets start simulating the transmission
	numSymbol = size(transmit_Mod, 2); % Number of transmitted symbols per stream per frame
	
    % Generate the channels. Assume channel to be independently fading
    % across symbols
    h = zeros(3, numSymbol); % The actual channels
    for i_transmission = 1 : 3
        h(i_transmission, :) = crandn(mu_h(i_transmission), sigma2_h(i_transmission), 1, numSymbol);
    end
    
    % Generate the received symbols
    y = zeros(2, numSymbol);
    y(1, :) = h(1, :) .* transmit_Mod(1, :);
    y(2, :) = sum(h(2 : 3, :) .* transmit_Mod(2 : 3, :), 1);
    y = y + crandn(0, sigma2_v, 2, numSymbol); % The received signals
    
   
    y = complex(real(y), imag(y)); % make sure y is complex, used for MAP
	h = complex(real(h), imag(h)); % make sure channel is complex, used for MAP
    
    %iterative receiver
    LextC = zeros(1, nldpc);
		
    error_perFrame_perMS = zeros(iter_max, 1);
    for iter = 1 : iter_max

        LextDemodulation = MAP_demod(y, h, bit_mat_anti, LextC, sym_mod_mat, sigma2_v);

        LextDemo_deinv = randdeintrlv(LextDemodulation, 0); %de-interleave
        [LLR_output_tmp, errors] = MpDecode(LextDemo_deinv, H_rows, H_cols, max_iterations, decoder_type, 1, 1, data); %ldpc decoder
        error_perFrame_perMS(iter) = errors(end); %count errors in each iteration
        LLR_output = LLR_output_tmp(end, :);
        LextDecoder = LLR_output - LextDemo_deinv; %calculate extrinic information from channel decoder
        LextC = randintrlv(LextDecoder, 0); %interleave
    end

    error_all(i) = error_perFrame_perMS(end);  

    %count errors
    if sum(error_all) > max_bit_error
        break;
    end
end
fprintf('\n');
codedBER = sum(error_all) / (bits_per_frame * i);

