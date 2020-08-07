function [am, bm, cm, dm] = Get_abcd_m( X1_next, X2_next, X1_current, X2_current, Hm, Gm2p, Gm1q, SCSI )
%GET_ABCD_M Summary of this function goes here
%   Detailed explanation goes here

if (nargin < nargin('Get_abcd_m'))
    SCSI = 0;
end

am = Get_ab_m(X1_next, X1_current, Hm, Gm2p, SCSI);
bm = Get_ab_m(X2_next, X2_current, Hm, Gm1q, SCSI);

cm = Get_cd_m(X1_current, Hm, Gm2p, SCSI);
dm = Get_cd_m(X2_current, Hm, Gm1q, SCSI);



end


function [ ab ] = Get_ab_m(X_next, X_current, Hm, Gm, SCSI)


DLU = size(X_next{1,1},2);
ULU = length(X_next{1,3});

temp = 0;

for iDLU = 1:1:DLU
    
    if (SCSI)
         temp = temp + real( (X_current{1,1}(:,iDLU))'*diag(diag(Hm*Hm'))*(X_next{1,1}(:,iDLU)) );
    else
    
        temp = temp + real( (X_current{1,1}(:,iDLU))'*Hm*Hm'*(X_next{1,1}(:,iDLU)) );
    end
    
end

if (SCSI)
    temp = temp + real( trace( (X_current{1,2})'*diag(diag(Hm*Hm'))*X_next{1,2} ) );
else
    temp = temp + real( trace( (X_current{1,2})'*Hm*Hm'*X_next{1,2} ) );
end

for iULU = 1:1:ULU
    
    temp = temp + X_current{1,3}(iULU)*X_next{1,3}(iULU)*(norm(Gm(iULU,:)))^2;
    
end

ab = temp;

end


function [ cd ] = Get_cd_m(X_current, Hm, Gm, SCSI)


DLU = size(X_current{1,1},2);
ULU = length(X_current{1,3});

temp = 0;

for iDLU = 1:1:DLU
    
    if (SCSI)
    
        temp = temp + (real((X_current{1,1}(:,iDLU))'*diag(diag(Hm*Hm'))*X_current{1,1}(:,iDLU) ));
    
    else
        temp = temp + (norm( Hm'*X_current{1,1}(:,iDLU) ))^2;
    end
    
end

if (SCSI)
    temp = temp + real(trace((X_current{1,2})'*diag(diag(Hm*Hm'))*X_current{1,2}) );
else
    temp = temp + (norm(Hm'*X_current{1,2},'fro') )^2;
end

for iULU = 1:1:ULU
    
    temp = temp + (X_current{1,3}(iULU))^2*(norm(Gm(iULU,:)))^2;
    
end

cd = temp;

end

