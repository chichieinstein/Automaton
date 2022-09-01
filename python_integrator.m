
function [iq, sample_rate] = python_integrator(label)
label = convertCharsToStrings(label);
rng(1);
RandHelper.setDefaultStream(1);
centerFreq_Hz       = 2.45e9; 
%% Source 1: WiFi

wlanCenterFreq_Hz = 2412e6; % Channel 1
nPacket = 6;
idleTime = [];
scramblerInitialization = [];
psduLength = 19;
message = [];

trafficType = atomic.Traffic('customArray','arrivalArray',[2 4 8 16 20 25 32 38 42]*1e-3); 

Wlan_HE_80211g = atomic.WlanNonHT80211g('wifi', trafficType, wlanCenterFreq_Hz, 1, nPacket, idleTime, scramblerInitialization, psduLength, message);


%% Source 2: BLE 

txPower_dbm         = 1;
bleCenterFreq_Hz    = 2426e6; 
message             = randi([0 1],640,1);
transmissionPerSec  = 100; 
trafficType         = atomic.Traffic('periodic', 'transmissionPerSec', transmissionPerSec); 
bleSymbolRate       = 1e6; 

ble_adv_channel_2426 = atomic.Bluetooth('ble_1',trafficType, bleCenterFreq_Hz, bleSymbolRate, txPower_dbm, message);

%% Source 3: LTE 

trafficType     = atomic.Traffic('periodic','transmissionPerSec',100); 
centerFreq_Hz_lte = 2.45e9 ; 
nPacket         = 0.1;
idleTime        = [];
TotSubframes    = 3.3;
message         = [];
NDLRB           = 6;

lte = atomic.LTE_DL_FDD('lte', trafficType, centerFreq_Hz_lte, 1, nPacket, idleTime, TotSubframes, message, NDLRB);

%% Source 4: DSSS

trafficType         = atomic.Traffic('periodic','transmissionPerSec',100); 
modOrder            = 2;
transmissionRate_Hz = 0.1e6;
samplesPerSymbol    = 2^7;
transmissionTotTime = 50e-3;
% freqSeparation_occupiedBandwidthFrac = 0.1;
spreadType          = 2;
samplesPerChip      = 10;
centerFreq_Hz_dsss  = 2.45e9;

ds3 =  atomic.Ds3('dsss', trafficType, centerFreq_Hz_dsss, transmissionRate_Hz, 10, samplesPerSymbol, modOrder, transmissionTotTime, spreadType,samplesPerChip);

%% Source 5: QAM 

trafficType = atomic.Traffic('periodic','transmissionPerSec',100); 
centerFreq_Hz_qam   = centerFreq_Hz ; 
transmissionRate_Hz = 5e6; 
samplesPerSymbol    = 2;
modOrder            = 16;
beta                = 0.4;
span                = 10;
transmissionTotTime = 5e-3;
txPowerdB           = 1;

qam                 = atomic.Qam('qam1',trafficType, centerFreq_Hz_qam, transmissionRate_Hz, txPowerdB, samplesPerSymbol, modOrder,beta, span,transmissionTotTime);

%% Source 6: PAM 

trafficType = atomic.Traffic('periodic','transmissionPerSec',100); 
centerFreq_Hz_pam   = 2.45e9; 
transmissionRate_Hz = 1e6; 
samplesPerSymbol    = 2;
modOrder            = 4;
beta                = 0.4;
span                = 10;
transmissionTotTime = 3.3e-3;
bandwidth_Hz_pam    = transmissionRate_Hz*(1+beta);
txPowerdB           = 1;

pam                 = atomic.Pam('pam1', trafficType, bandwidth_Hz_pam, centerFreq_Hz_pam, transmissionRate_Hz, txPowerdB, modOrder, samplesPerSymbol, beta, span, transmissionTotTime);

%% Source 7: FSK 

trafficType = atomic.Traffic('periodic','transmissionPerSec',100); 
modOrder            = 2;
transmissionRate_Hz = 2e6;
samplesPerSymbol    = 2;
transmissionTotTime = 2e-3;
freqDev             = 0.5;
txpowerdB           = 1;
centerFreq_Hz_fsk   = 2.45e9;

fsk = atomic.Fsk('fsk', trafficType, centerFreq_Hz_fsk, transmissionRate_Hz, samplesPerSymbol, modOrder, txpowerdB, freqDev,transmissionTotTime);

%% Source 8: GFSK

trafficType = atomic.Traffic('periodic','transmissionPerSec',100); 
modOrder            = 4;
transmissionRate_Hz = 5e6;
samplesPerSymbol    = 2;
transmissionTotTime = 2e-3;
txpowerdB           = 1;
bandwidth_Hz_gfsk   = transmissionRate_Hz*2;
centerFreq_Hz_gfsk = 2.45e9-35e6;

gfsk = atomic.Gfsk('gfsk_1', trafficType, bandwidth_Hz_gfsk, centerFreq_Hz_gfsk, transmissionRate_Hz, txpowerdB, samplesPerSymbol, modOrder,transmissionTotTime);

%% Source 9: PSK 

trafficType = atomic.Traffic('periodic','transmissionPerSec',100); 
modOrder            = 4;
transmissionRate_Hz = 2e6; 
samplesPerSymbol    = 4;
transmissionTotTime = 1e-3;
beta                = 0.4;
span                = 10;
txpowerdB           = 1;
centerFreq_Hz_psk   = 2.45e9;

psk                 = atomic.Psk('psk_1', trafficType, centerFreq_Hz_psk, transmissionRate_Hz, txpowerdB, samplesPerSymbol, modOrder,beta,span,transmissionTotTime);


%% Source 10: OFDM 

trafficType = atomic.Traffic('periodic','transmissionPerSec',100); 
modOrder                = 4;
transmissionRate_Hz     = 8e6;
transmissionTotalTime   = 4e-3;
Nsc                     = 64;
tx_powerdb              = 1;
centerFreq_Hz_ofdm      = 2.45e9;

ofdm                    =  atomic.Ofdm("ofdm_1", trafficType, centerFreq_Hz_ofdm, transmissionRate_Hz,transmissionTotalTime,tx_powerdb, Nsc, modOrder);

%% Defining the rx in the end so no resample occurs

tot_time            = 5e-3; % same as the one you specified for the signal

% Defining rx sample rate same as the signal's sample rate so no resample occurs
% rxSampleRate_Hz     = ds3.transmissionRate_Hz; % signal name you defined above_transmissionRate_Hz, 
                                               % eg: qam.transmissionRate_Hz

%% Choose the signal you want (keep the rest commented)

rxLocation          = [0 0 0];
channel = "RICIANNORM";
sourceName = 'source1';
sourceLoc = (2*rand(1,3)-1)*10;
sourceFreqOffset = 0; % realistic CFO is 1 ppm of Fc
imperfectionCfg = atomic.RFImperfections(sourceFreqOffset);

switch label
    case "wlan"
        rxSampleRate_Hz = Wlan_HE_80211g.transmissionRate_Hz;
        source = atomic.Source(sourceName, 'custom label', rxSampleRate_Hz, sourceLoc, channel, imperfectionCfg);
        source.addSignal(Wlan_HE_80211g);
        rx = atomic.Rx('rx1',rxSampleRate_Hz, wlanCenterFreq_Hz, rxLocation);
        var = Wlan_HE_80211g;
    case "ofdm"
        rxSampleRate_Hz = ofdm.transmissionRate_Hz;
        source = atomic.Source(sourceName, 'custom label', rxSampleRate_Hz, sourceLoc, channel, imperfectionCfg);
        source.addSignal(ofdm);
        rx = atomic.Rx('rx1',rxSampleRate_Hz, centerFreq_Hz_ofdm, rxLocation);
        var = ofdm;
    case "psk"
        rxSampleRate_Hz = psk.transmissionRate_Hz;
        source = atomic.Source(sourceName, 'custom label', rxSampleRate_Hz, sourceLoc, channel, imperfectionCfg);
        source.addSignal(psk);
        rx = atomic.Rx('rx1',rxSampleRate_Hz, centerFreq_Hz_psk, rxLocation);
        var = psk;
    case "gfsk"
        rxSampleRate_Hz = gfsk.transmissionRate_Hz;
        source = atomic.Source(sourceName, 'custom label', rxSampleRate_Hz, sourceLoc, channel, imperfectionCfg);
        source.addSignal(gfsk);
        rx = atomic.Rx('rx1',rxSampleRate_Hz, centerFreq_Hz_gfsk, rxLocation);
        var = gfsk;
    case "fsk"
        rxSampleRate_Hz = fsk.transmissionRate_Hz;
        source = atomic.Source(sourceName, 'custom label', rxSampleRate_Hz, sourceLoc, channel, imperfectionCfg);
        source.addSignal(fsk);
        rx = atomic.Rx('rx1',rxSampleRate_Hz, centerFreq_Hz_fsk, rxLocation);
        var = fsk;
    case "pam"
        rxSampleRate_Hz = pam.transmissionRate_Hz;
        source = atomic.Source(sourceName, 'custom label', rxSampleRate_Hz, sourceLoc, channel, imperfectionCfg);
        source.addSignal(pam);
        rx = atomic.Rx('rx1',rxSampleRate_Hz, centerFreq_Hz_pam, rxLocation);
        var = pam;
    case "qam"
        rxSampleRate_Hz = qam.transmissionRate_Hz;
        source = atomic.Source(sourceName, 'custom label', rxSampleRate_Hz, sourceLoc, channel, imperfectionCfg);
        source.addSignal(qam);
        rx = atomic.Rx('rx1',rxSampleRate_Hz, centerFreq_Hz_qam, rxLocation);
        var = qam;
    case "ds3"
        rxSampleRate_Hz =ds3.transmissionRate_Hz;
        source = atomic.Source(sourceName, 'custom label', rxSampleRate_Hz, sourceLoc, channel, imperfectionCfg);
        source.addSignal(ds3);
        rx = atomic.Rx('rx1',rxSampleRate_Hz, centerFreq_Hz_dsss, rxLocation);
        var = ds3;
    case "lte"
        rxSampleRate_Hz = lte.transmissionRate_Hz;
        source = atomic.Source(sourceName, 'custom label', rxSampleRate_Hz, sourceLoc, channel, imperfectionCfg);
        source.addSignal(lte);
        rx = atomic.Rx('rx1',rxSampleRate_Hz, centerFreq_Hz_lte, rxLocation);
        var = lte;
    case "ble"
        rxSampleRate_Hz = ble_adv_channel_2426.transmissionRate_Hz;
        source = atomic.Source(sourceName, 'custom label', rxSampleRate_Hz, sourceLoc, channel, imperfectionCfg);
        source.addSignal(ble_adv_channel_2426);
        rx = atomic.Rx('rx1',rxSampleRate_Hz, bleCenterFreq_Hz, rxLocation);
        var = ble_adv_channel_2426;
end

bandwidth_Hz        = 1*rxSampleRate_Hz;

%% Thermal nois
sourceThermal        = atomic.Source('Noise', 'custom label',rxSampleRate_Hz, [0 0 0]);
thermalTemp_K        = 1;

wb_thermal_noise_hdl = atomic.WidebandThermalWgn(...
                                                    "thermal_1",...
                                                    bandwidth_Hz, ...
                                                    centerFreq_Hz, ...
                                                    rxSampleRate_Hz, ...
                                                    thermalTemp_K);

sourceThermal.addSignal(wb_thermal_noise_hdl);

%% Running the dataset generator

folder          = "/home/srajagopal/Desktop/scisrs_dataset/" + label;
filename_base   = 'data_' + string(convertTo(datetime, "yyyymmdd")) + '_' + string(rxSampleRate_Hz);
sigsPerEsNo         = 1;

EsNoRange           = 10;

% Give the name of the signal you desire as the first argument 
txPowdBmVec = wb_thermal_noise_hdl.getTxPowVecFromEsNo(var, EsNoRange);
filename_base_iter = filename_base + "_" + string(EsNoRange(1)) +"dB";
sigGen = VirtualSignalEngine();
sigGen.addRxObj(rx);
sigGen.addSource(source);
sigGen.addSource(sourceThermal);
tic
sigGen.sourceArray(1).signalArray(1).setTxPowdBm(txPowdBmVec(1));
sigGen.generateMultipleDataSets(0, tot_time, sigsPerEsNo, folder, filename_base_iter); %10 files per EsNo
toc

dataPointIdx = 1;
EsNo = EsNoRange(1);

[readSamplesIQ, readMetadataStruct] = VirtualSignalEngine.readDataFiles(folder,...
            string(filename_base)+"_"+string(EsNo) +"dB"+"_"+string(dataPointIdx));
sample_rate = rxSampleRate_Hz;
iq = readSamplesIQ;

end