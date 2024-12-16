function [lg2,fou2] = gen_doubleW(lg,fou,nvein)

% Generates lg and fou structures, as those generated with file2body an file2fou functions, in a format
% to be read by tucan launcher with double wall. For that, two aditional bodies full of zeros are included
% at the end of the arrays with the length of LA and MV arrays

lg2 = lg;
fou2 = fou;

% MV body index
iMV = nvein+2;

% number of points of LA and MV
nLA = lg.ie(1);
nMV = lg.ie(nvein+2)-lg.ie(nvein+1);

% vector full of zeros to be added at the end of the arrays
vzeros = zeros(nLA+nMV,1);

% update the fields of lg struct
lg2.nreal = lg.nreal+nLA+nMV;
lg2.xyz = [lg.xyz; [vzeros vzeros vzeros] ];
lg2.vol = [lg.vol; vzeros];
lg2.normvec = [lg.normvec; [vzeros vzeros vzeros] ];
lg2.xyzc = [lg.xyzc; lg.xyzc(1,:); lg.xyzc(iMV,:)];
lg2.ib = [lg.ib; lg.ie(end)+1; lg.ie(end)+1+nLA];
lg2.ie = [lg.ie; lg.ie(end)+nLA; lg.ie(end)+nLA+nMV];

ndim = lg.ndim;
nspi = 0.5*ndim*(ndim+1);
lg2.spi = [lg.spi; lg.spi(1:nspi); lg.spi(((iMV-1)*nspi+1):iMV*nspi)];

% update the field of fou struct
fou2.nreal = fou.nreal+nLA+nMV;

fou2.xyz_r = zeros(fou.nreal+nLA+nMV,3,fou.nm1);
fou2.xyz_i = zeros(fou.nreal+nLA+nMV,3,fou.nm1);
fou2.xyz_r(1:fou.nreal,:,:) = fou.xyz_r;
fou2.xyz_i(1:fou.nreal,:,:) = fou.xyz_i;

fou2.vecprod_r = zeros(fou.nreal+nLA+nMV,3,fou.nm2);
fou2.vecprod_i = zeros(fou.nreal+nLA+nMV,3,fou.nm2);
fou2.vecprod_r(1:fou.nreal,:,:) = fou.vecprod_r;
fou2.vecprod_i(1:fou.nreal,:,:) = fou.vecprod_i;

end
