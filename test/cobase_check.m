bas_out;

function [S] = bas_submat(Mat, Cob)
	[n_m, d_m] = size(Mat);
	[n_c, d_c] = size(Cob);
	S = zeros(d_c, d_m);
	for i = 1:d_c
		S(i,:) = Mat(Cob(i), :);
	end;
endfunction;

function [nRanks] = bas_check_ranks(Mat, Cobs)
	nRanks = 0;
	[n, d] = size(Cobs);
	for i = 1:n
		S = bas_submat(Mat, Cobs(i,:));
		r = rank(S);
		disp(r);
		if ( r ~= d ) nRanks += 1; end;
	end;
endfunction;

disp(bas_check_ranks(M, C));
