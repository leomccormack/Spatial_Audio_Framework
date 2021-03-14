clear all, dbstop if error, close all %#ok

% afSTFT

 % 0.070644 seconds. to do this with the C interface
nCHin = 6;
nCHout = 13;
hopsize = 128;
blocksize = 2048;


mh = mexhost;
[freqVector, procDelay] = feval(mh,'saf_afSTFT', nCHin, nCHout, hopsize, 1, 0, 48e3);
tic
for i=1:round(48000/blocksize * 2)
    insig = randn(blocksize, nCHin);  
    outsig = feval(mh,'saf_afSTFT', insig);  
     outsig = repmat(outsig(:,1,:), [1 nCHout 1 ]);
     testsig = feval(mh,'saf_afSTFT', outsig);  
end
toc
feval(mh,'saf_afSTFT')

 
[freqVector, procDelay] = saf_afSTFT(nCHin, nCHout, hopsize, 1, 0, 48e3);
tic
for i=1:round(48000/blocksize * 2)
    insig = randn(blocksize, nCHin);  
    outsig = saf_afSTFT(insig);  
    outsig = repmat(outsig(:,1,:), [1 nCHout 1 ]);
    testsig = saf_afSTFT(outsig);  
end
toc
saf_afSTFT();
 
 

% [freqVector, procDelay] = safmex_afSTFT(nCHin, nCHout, hopsize, blocksize, 1, 1, 48e3); 
% dataTD_ref = randn(blocksize, nCHin); 
% dataFD = safmex_afSTFT(dataTD_ref.');
% dataFD = repmat(dataFD(:,1,:), [1 nCHout 1]); % copy 1st input channel to all output channels
% dataTD = safmex_afSTFT(dataFD);
% dataTD_ref = dataTD_ref(1:end-procDelay,1);
% dataTD = dataTD(1,procDelay+1:end).';
% nTests = nTests+1;
% if max(abs(dataTD(:,1)-dataTD_ref(:,1)))<0.01, nPass=nPass+1; else, nFail=nFail+1; end 
% safmex_afSTFT();



% 
% 
% %% TEST CONFIG
% nTests = 0; nPass = 0; nFail = 0;
% tol = 1e-5; % FLT_EPSILON
%  
% %% TESTS
% % safmex_tracker3d 
% nTests = nTests+1; 
% safmex_tracker3d_tests  
% 
% % safmex_latticeDecorrelator
% nCH = 6; 
% hopsize = 128;
% blocksize = 2048*24;
% [freqVector, procDelay] = safmex_afSTFT(nCH, nCH, hopsize, blocksize, 1, 0, 48e3); 
% orders = [6 6 4 2].';
% fixedDelays = [ 8 6 4 2 2].';
% freqCutoffs = [ 300 800 2e3, 4e3].';
% timeslots = blocksize/hopsize;
% safmex_latticeDecorrelator(nCH, orders, freqCutoffs, fixedDelays, freqVector, timeslots);
% dataTD_ref = randn(blocksize, nCH); 
% dataFD = safmex_afSTFT(dataTD_ref.');
% dataFD = safmex_latticeDecorrelator(dataFD);
% dataTD = safmex_afSTFT(dataFD);
% dataTD_ref = dataTD_ref(1:end-procDelay,1);
% dataTD = dataTD(1,procDelay+1:end).';
% icc = (dataTD_ref.' * dataTD) / sqrt((dataTD_ref.' * dataTD_ref) * (dataTD.' * dataTD));
% nTests = nTests+1;
% if icc<0.05, nPass=nPass+1; else, nFail=nFail+1; end 
% safmex_afSTFT();
% safmex_latticeDecorrelator();
% 
% % safmex_faf_IIRFilterbank 
% order = 3;
% lSig = 2048*10;
% fs = 48e3;
% cutoffFreqs = [125 250 500 1000 2000 4000]*2/sqrt(2);
% safmex_faf_IIRFilterbank(order, cutoffFreqs.', lSig, fs); 
% data_ref = zeros(lSig, 1);  data_ref(1,1) = 1;
% data_bands = safmex_faf_IIRFilterbank(data_ref).';
% figure;  
% f = linspace(1, fs / 2, lSig / 2 + 1 );
% for idx = 1:size(data_bands, 2)
%     vTf = fft(data_bands(:, idx)) ./ fft(data_ref);
%     vAbsTf = 20 * log10(abs(vTf(1:end / 2 + 1)));
%     semilogx(f, vAbsTf), hold on
% end
% vYsum = sum(data_bands, 2);
% vTfSum = fft(vYsum) ./ fft(data_ref);
% vAbsTfSum = 20 * log10(abs(vTfSum(1:end / 2 + 1)));
% semilogx(f, vAbsTfSum, 'k'), hold on
% xlabel('Frequency [Hz]'), ylabel('Magnitude [dB]'), title('sum(bands) should be ~0dB')
% ylim([-60 10]), xlim([20 20e3]), grid on, hold off
% nTests = nTests+1;
% if max(abs(vAbsTfSum(:)))<0.1, nPass=nPass+1; else, nFail=nFail+1; end 
% safmex_faf_IIRFilterbank();
% 
% % safmex_afSTFT
% nCHin = 6;
% nCHout = 13;
% hopsize = 128;
% blocksize = 2048*24;
% [freqVector, procDelay] = safmex_afSTFT(nCHin, nCHout, hopsize, blocksize, 1, 1, 48e3); 
% dataTD_ref = randn(blocksize, nCHin); 
% dataFD = safmex_afSTFT(dataTD_ref.');
% dataFD = repmat(dataFD(:,1,:), [1 nCHout 1]); % copy 1st input channel to all output channels
% dataTD = safmex_afSTFT(dataFD);
% dataTD_ref = dataTD_ref(1:end-procDelay,1);
% dataTD = dataTD(1,procDelay+1:end).';
% nTests = nTests+1;
% if max(abs(dataTD(:,1)-dataTD_ref(:,1)))<0.01, nPass=nPass+1; else, nFail=nFail+1; end 
% safmex_afSTFT();
% 
% % safmex_qmf
% nCHin = 10;
% nCHout = 4;
% hopsize = 128;
% blocksize = 2048*20;
% [freqVector2, procDelay] = safmex_qmf(nCHin, nCHout, hopsize, blocksize, 1, 0, 48e3); 
% dataTD_ref = randn(blocksize, nCHin); 
% dataFD = safmex_qmf(dataTD_ref.');
% dataFD = repmat(dataFD(:,1,:), [1 nCHout 1]); % copy 1st input channel to all output channels
% dataTD = safmex_qmf(dataFD);
% dataTD_ref = dataTD_ref(1:end-procDelay,1);
% dataTD = dataTD(1,procDelay+1:end).';
% nTests = nTests+1;
% if max(abs(dataTD(:,1)-dataTD_ref(:,1)))<0.01, nPass=nPass+1; else, nFail=nFail+1; end 
% safmex_qmf();
% 
% % safmex_generateVBAPgainTable3D
% [~,ls_dirs] = getTdesign(10);
% ls_dirs = ls_dirs*180/pi;
% aziElevRes = [5 5];
% gtable = safmex_generateVBAPgainTable3D(ls_dirs, aziElevRes(1), aziElevRes(2), 1, 0, 0);
% gtable_ref = getGainTable(ls_dirs,aziElevRes,0,'vbap');
% nTests = nTests+1;
% if max(abs(gtable_ref(:)-gtable(:)))<tol, nPass=nPass+1; else, nFail=nFail+1; end 
% 
% % safmex_getSHreal
% order = 5;
% dirs = randn(1000,2); 
% Y = safmex_getSHreal(order,dirs); 
% Y_ref = getSH(order, dirs, 'real').'; 
% nTests = nTests+1;
% if max(abs(Y(:)-Y_ref(:)))<tol, nPass=nPass+1; else, nFail=nFail+1; end 
% 
% % safmex_getSHcomplex
% order = 4;
% dirs = randn(400,2); 
% Y = safmex_getSHcomplex(order,dirs); %Y = 
% Y_ref = getSH(order, dirs, 'complex').';
% nTests = nTests+1;
% if max(abs(Y(:)-Y_ref(:)))<tol, nPass=nPass+1; else, nFail=nFail+1; end 
% 
% % Output
% display(['Number of tests: ' num2str(nTests)]);
% display(['Number of failures: ' num2str(nFail)]);
