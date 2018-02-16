function [ pak , fo_sp ] = fFD_preamble_detect_final(Fs,samples,short_preamble,threshold)

s=length(samples);
proc=zeros(s,1);
b=fft([short_preamble;zeros(352,1)]);
tocopy = 8000;
toskip = 5;
pak=zeros(tocopy,1);
i=1;
j=1;
K=floor(s/353);
x=zeros(K,1);
conv_prev=zeros(512,1);
skip_i=0;
fo_sp=zeros(1,1);
while i<=K-toskip
   batch=[samples((i-1)*353+1:i*353);zeros(159,1)]; 
   a=fft(batch);
   c=a.*b;
   d=ifft(c);
   d(1:159)=d(1:159)+conv_prev(end-158:end);
   proc((i-1)*353+1:i*353)=abs(d(1:353));
   
   [x,idx]=max(abs(d(1:353))>threshold);
   if ~isempty(x) && x~=0 && i>skip_i
      skip_i=i+toskip;
      pak(:,j)=samples((i-1)*353+idx:(i-1)*353+tocopy+idx-1);

      % Frequency offset correction using the Short Preamble
      [fo_sp(1,j)] = fFD_frequency_offset_estimation(Fs,pak(:,j),...
                    'short_preamble');
      [pak(:,j)] = fFD_frequency_offset(Fs,-fo_sp(1,j),pak(:,j));
      j=j+1;
   end
   i=i+1; 
   conv_prev=d;
end

end