function [ packet, lptot, fo_lp, num_pkt ] = fFD_symbolalign2(Fs,pak,tx_preamble)
num_pkt=size(pak,2);
lp1=zeros(64,1);
lp2=zeros(64,1);
lptot = zeros(160,1);
tocopy = 8000;
packet=zeros(tocopy,num_pkt);

cor = zeros(512,num_pkt);
fo_lp = zeros(1,1);
for jj=1:num_pkt

preamble_area=pak(1:449,jj);
bb=fft(tx_preamble,512);
aa=fft(preamble_area,512);
cc=aa.*bb;
corr_out=abs(ifft(cc));
cor(:,jj) = corr_out;

[peak1,indx1]=max(corr_out);

corr_out(indx1)=0;

[peak2,indx2]=max(corr_out);

corr_out(indx2)=0;

[peak3,indx3]=max(corr_out);

if ~(abs(indx1-indx2)==64 && (indx1-indx3==64 || indx2-indx3==64))
%     display('right peaks not found!');
    %continue
end

if indx1 > indx2
    start_frame_index=indx1+1;
    found = 1;
    %start_frame_index=511-indx1+64;
%     display('index1 found!')
else
    start_frame_index=indx2+1;    
    found = 1;
    %start_frnum_pktame_index=511-indx2+64;
%     display('index1 found!')
end

% longpre
if(start_frame_index>160)
    packet(:,jj)=[pak(start_frame_index:end,jj);zeros(tocopy-length(pak(start_frame_index:end,jj)),1)];
    lp1(:,jj)=pak(start_frame_index-64:start_frame_index-1,jj);
    lp2(:,jj)=pak(start_frame_index-64-64:start_frame_index-64-1,jj);
    lptot(:,jj)=pak(start_frame_index-64-64-32:start_frame_index-1,jj);

    % Frequency offset correction using the Long Preamble
    [fo_lp(1,jj)] = fFD_frequency_offset_estimation(Fs,lptot(:,jj),'long_preamble');
    [lptot(:,jj)] = fFD_frequency_offset(Fs,-fo_lp(1,jj),lptot(:,jj));
    [packet(:,jj)] = fFD_frequency_offset(Fs,-fo_lp(1,jj),packet(:,jj));
end
end

end
