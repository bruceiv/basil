img_mat;
cob_imgs;

[n, d] = size(M);

function [aMat] = aug_mat(Mat)
	[n_m, d_m] = size(Mat);
	aMat = [ Mat ; 1 zeros(1,d_m-1) ];
endfunction;

B = inv(aug_mat(M(s,:)));

[vn, vd] = size(vs);

for i = 1:vn
	v = vs(i,:);
	T = B * aug_mat(M(v,:));
	Mt = (T * M')';
	if ( all(M == Mt) )
		disp(["PASS     " gs(i,:)]);
	else
		disp(["    FAIL " gs(i,:)]);
		disp([ repmat("         ",n,1) num2str(M == Mt) repmat(" | ",n,1) num2str(M) repmat(" | ",n,1) num2str(Mt) ]);
	end;
end;
