clear all;close all;
% rng(0);
% CPM parameter
h = 0.5; % modulation number
L = 2;   % related length
nn = 50; % No. of modulate code
Tp = 5;  % LFM time width
Ts = Tp/nn;% symbol time width
sym_rate = 1/Ts;
% LFM parameter
B = 200; % Band width MHz
K = B/Tp;  % FM slope MHz/mus B/Tp
fc = 200; % carry frequence MHz 
fd = 0.5;
phid = 0.1 * pi;
fs = 640; % sample frequence MHz 
oversample = fs*Ts;
%fs = 1600;
Ne = 5000;

Co_rec = [];
D1_rec = [];
D2_rec = [];

for EbN0 = 0:10

ACo = 0;
AD1 = 0;
AD2 = 0;
erro_cnt_Co = 0; % statistic erro rate
erro_cnt_D1 = 0; % statistic erro rate
erro_cnt_D2 = 0; % statistic erro rate

for ll = 1:Ne
data_pre = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
data = randi([0 1], [1 nn-15]);
CPM_BB = CPMmod([data_pre, data], oversample);
% figure
% plot(angle(CPM_BB))
% figure
% plot(angle(CPM_BB(16:16:end)))
% hold on
% figure
% plot((10*log10(abs(fft(CPM_BB(1:16:end))))))
% hold on
t = -0.5*Tp : 1/fs: 0.5*Tp - 1/fs;
t = t + 1/fs/2;

CPM_LFM = CPM_BB .* exp(1i*K*pi*t.^2) .* exp(1i*2*pi*fc*t) .* exp(1i*phid);
% CPM_LFM = CPM_BB;
% figure
% plot(abs((fft(CPM_LFM))))
% hold on;
CPM_LFM_noise = awgn(CPM_LFM, -10*log10(64) + EbN0, 'measured');

% 640MHz 90 110 290 310
BP = [0.000166953525860736,-0.000709394692822494,-0.00123272880128526,0.000389568141776851,-7.93908734925980e-05,0.00201400316864084,-0.000693442506475319,0.000606135454402551,-0.00320390518367515,0.000734938761785564,-0.00105362130244274,0.00407432957973307,8.92631027491000e-05,0.000924995775111281,-0.00400124990996845,-0.00233873293791828,0.000120264204737338,0.00256995425101018,0.00619792158218013,-0.00194245677490267,9.26554636924372e-05,-0.0111579741746940,0.00365902103596419,-0.00312353218053580,0.0159073275425811,-0.00364103718212918,0.00504633440967609,-0.0185673731733380,-0.000136658472021936,-0.00426207697051575,0.0172772538514604,0.00947315623028422,-0.000184567415454268,-0.0109421093376432,-0.0252962452031788,0.00779423683768114,-0.000139034382431863,0.0474168990552993,-0.0158877089214669,0.0140749417827033,-0.0753680426873037,0.0184918870429166,-0.0277619144713638,0.112828785028356,0.000170328705549544,0.0377707932992514,-0.207029038134832,-0.202625979602336,0.625225670590477,-0.202625979602336,-0.207029038134832,0.0377707932992514,0.000170328705549544,0.112828785028356,-0.0277619144713638,0.0184918870429166,-0.0753680426873037,0.0140749417827033,-0.0158877089214669,0.0474168990552993,-0.000139034382431863,0.00779423683768114,-0.0252962452031788,-0.0109421093376432,-0.000184567415454268,0.00947315623028422,0.0172772538514604,-0.00426207697051575,-0.000136658472021936,-0.0185673731733380,0.00504633440967609,-0.00364103718212918,0.0159073275425811,-0.00312353218053580,0.00365902103596419,-0.0111579741746940,9.26554636924372e-05,-0.00194245677490267,0.00619792158218013,0.00256995425101018,0.000120264204737338,-0.00233873293791828,-0.00400124990996845,0.000924995775111281,8.92631027491000e-05,0.00407432957973307,-0.00105362130244274,0.000734938761785564,-0.00320390518367515,0.000606135454402551,-0.000693442506475319,0.00201400316864084,-7.93908734925980e-05,0.000389568141776851,-0.00123272880128526,-0.000709394692822494,0.000166953525860736];
% 640MHz 90 110
LP = [-0.000304477932586780,0.00102468017652464,0.000780222316046232,0.000196446608140970,-0.000660421418668351,-0.00107247326501243,-0.000483410134579629,0.000790121227889237,0.00161919883369055,0.00101907779873649,-0.000780329818227999,-0.00225160180996632,-0.00181591831693872,0.000549685538050848,0.00291562973304411,0.00290842422179811,7.37206174068678e-06,-0.00352422810673698,-0.00431133033600945,-0.00100627640660450,0.00395557431229063,0.00601226653687585,0.00257065714222627,-0.00405219447425498,-0.00797019412674347,-0.00483492316513782,0.00361284211665321,0.0101133597346884,0.00795466416012556,-0.00237932839477009,-0.0123445931136466,-0.0121446190232754,-1.17644278880950e-05,0.0145451857816449,0.0177716919364819,0.00414817157611159,-0.0165853807141233,-0.0256129514393360,-0.0112035632931889,0.0183368219677982,0.0376945155969487,0.0242509155866588,-0.0196823884424215,-0.0611356108881986,-0.0555265271475656,0.0205304488837567,0.146540933939029,0.264450071755388,0.312513994148623,0.264450071755388,0.146540933939029,0.0205304488837567,-0.0555265271475656,-0.0611356108881986,-0.0196823884424215,0.0242509155866588,0.0376945155969487,0.0183368219677982,-0.0112035632931889,-0.0256129514393360,-0.0165853807141233,0.00414817157611159,0.0177716919364819,0.0145451857816449,-1.17644278880950e-05,-0.0121446190232754,-0.0123445931136466,-0.00237932839477009,0.00795466416012556,0.0101133597346884,0.00361284211665321,-0.00483492316513782,-0.00797019412674347,-0.00405219447425498,0.00257065714222627,0.00601226653687585,0.00395557431229063,-0.00100627640660450,-0.00431133033600945,-0.00352422810673698,7.37206174068678e-06,0.00290842422179811,0.00291562973304411,0.000549685538050848,-0.00181591831693872,-0.00225160180996632,-0.000780329818227999,0.00101907779873649,0.00161919883369055,0.000790121227889237,-0.000483410134579629,-0.00107247326501243,-0.000660421418668351,0.000196446608140970,0.000780222316046232,0.00102468017652464,-0.000304477932586780];

CPM_LFM_recv = conv(CPM_LFM_noise, BP);
CPM_LFM_recv = CPM_LFM_recv(48+1:end-48);

% plot(10*log10(abs((fft(CPM_LFM_recv)))))

CPM_recv = CPM_LFM_recv .* exp(-1i*K*pi*t.^2) .* exp(-1i*2*pi*(fc-fd)*t);
% CPM_recv = CPM_LFM_recv;
% figure
% plot(10*log10(abs((fft(CPM_recv)))))
% 640MHz 20 50
LP_base1 = [-9.72085947313688e-05,-4.06529358409810e-05,-1.92706267281707e-05,3.07331626878567e-05,0.000111159511298551,0.000217033647458863,0.000334807463032157,0.000442234611454024,0.000509982842655053,0.000505507124986668,0.000398783267474579,0.000169485465682687,-0.000185678523900577,-0.000646750347963097,-0.00116707310421438,-0.00167387344460007,-0.00207395653711537,-0.00226458781510777,-0.00214904487736866,-0.00165503145373915,-0.000753560872558373,0.000525106467162800,0.00208153056379269,0.00374630922014435,0.00529022407301724,0.00644637958023135,0.00694280842171609,0.00654282801915436,0.00508800890259948,0.00253862170635199,-0.000995314637410876,-0.00523815427766679,-0.00975265055446563,-0.0139655069552991,-0.0172148628042430,-0.0188161766424628,-0.0181389057992720,-0.0146876361500108,-0.00817523066798105,0.00141948185535040,0.0138200002949033,0.0284508547993113,0.0444702674437088,0.0608342564632763,0.0763872448444627,0.0899691719967881,0.100528110346876,0.107225111712627,0.109519659248906,0.107225111712627,0.100528110346876,0.0899691719967881,0.0763872448444627,0.0608342564632763,0.0444702674437088,0.0284508547993113,0.0138200002949033,0.00141948185535040,-0.00817523066798105,-0.0146876361500108,-0.0181389057992720,-0.0188161766424628,-0.0172148628042430,-0.0139655069552991,-0.00975265055446563,-0.00523815427766679,-0.000995314637410876,0.00253862170635199,0.00508800890259948,0.00654282801915436,0.00694280842171609,0.00644637958023135,0.00529022407301724,0.00374630922014435,0.00208153056379269,0.000525106467162800,-0.000753560872558373,-0.00165503145373915,-0.00214904487736866,-0.00226458781510777,-0.00207395653711537,-0.00167387344460007,-0.00116707310421438,-0.000646750347963097,-0.000185678523900577,0.000169485465682687,0.000398783267474579,0.000505507124986668,0.000509982842655053,0.000442234611454024,0.000334807463032157,0.000217033647458863,0.000111159511298551,3.07331626878567e-05,-1.92706267281707e-05,-4.06529358409810e-05,-9.72085947313688e-05];
CPM_recv = conv(CPM_recv, LP_base1);
CPM_recv = CPM_recv(48+1:end-48);
CPM_recv = CPM_recv(1:4:end);
% figure
% plot((10*log10(abs(fft(CPM_recv)))))
% 160MHz 10 16
LP_base2 = [4.36819628576457e-06,-0.000466221503953210,-0.000384243939491384,-0.000386666869176788,-0.000250635733261967,3.45947659446817e-05,0.000425341094721902,0.000823112583152469,0.00109251835769622,0.00109766663112481,0.000749880051849361,5.25395790953453e-05,-0.000873401826615632,-0.00179716610534241,-0.00242468901921108,-0.00247831995734417,-0.00179062733633171,-0.000385117810334547,0.00148496511257771,0.00336229382178751,0.00468059432117910,0.00491486938309890,0.00375080866702495,0.00122570895358806,-0.00220860482157948,-0.00573428412476804,-0.00833592147844058,-0.00905864444529901,-0.00730161401124888,-0.00305991282536211,0.00295339670466841,0.00937202143866384,0.0144348545004435,0.0164011655097522,0.0140383048172512,0.00706052460461645,-0.00360438651221782,-0.0158126071026934,-0.0265040719900253,-0.0322711849339532,-0.0300897557984063,-0.0180548029207904,0.00404990014886409,0.0344962487869823,0.0697552798215495,0.105034864314983,0.135118587339890,0.155335467890678,0.162458306665658,0.155335467890678,0.135118587339890,0.105034864314983,0.0697552798215495,0.0344962487869823,0.00404990014886409,-0.0180548029207904,-0.0300897557984063,-0.0322711849339532,-0.0265040719900253,-0.0158126071026934,-0.00360438651221782,0.00706052460461645,0.0140383048172512,0.0164011655097522,0.0144348545004435,0.00937202143866384,0.00295339670466841,-0.00305991282536211,-0.00730161401124888,-0.00905864444529901,-0.00833592147844058,-0.00573428412476804,-0.00220860482157948,0.00122570895358806,0.00375080866702495,0.00491486938309890,0.00468059432117910,0.00336229382178751,0.00148496511257771,-0.000385117810334547,-0.00179062733633171,-0.00247831995734417,-0.00242468901921108,-0.00179716610534241,-0.000873401826615632,5.25395790953453e-05,0.000749880051849361,0.00109766663112481,0.00109251835769622,0.000823112583152469,0.000425341094721902,3.45947659446817e-05,-0.000250635733261967,-0.000386666869176788,-0.000384243939491384,-0.000466221503953210,4.36819628576457e-06];
CPM_recv = conv(CPM_recv, LP_base2);
CPM_recv = CPM_recv(48+1:end-48);
CPM_recv = CPM_recv(1:4:end);
% figure
% plot((10*log10(abs(fft(CPM_recv)))))
% 40MHz 9.5 11
% 40MHz 5 7
% 4 6
% 5.5 7.5
% 4.5 6
% LP_base3 = [0.000330237662212286,0.000198295103017753,-0.000469504941879963,-3.07332198143413e-05,0.000530446929797254,0.000110049437588789,-0.000838089356596170,-1.59234941706467e-05,0.00109889241302976,-4.63448770326659e-05,-0.00147731989791978,0.000209215904988000,0.00189134412292593,-0.000428686650077742,-0.00238769385139514,0.000758811266874641,0.00294151183186525,-0.00120602897645457,-0.00356275247645879,0.00180683072389220,0.00423922975437062,-0.00258731638098111,-0.00496649404133663,0.00358652801599187,0.00573261054482327,-0.00484707120349307,-0.00652616319850337,0.00642411204097324,0.00733168780476571,-0.00838972906356833,-0.00813404310905357,0.0108448452085734,0.00891540962039221,-0.0139392187955169,-0.00965872083006844,0.0179123452201005,0.0103462462982944,-0.0231767317862166,-0.0109610400592562,0.0305139920646979,0.0114877299726374,-0.0416075114003120,-0.0119125466162078,0.0608515715420498,0.0122244874648869,-0.104396997365306,-0.0124149646578933,0.317737757575639,0.512478997233137,0.317737757575639,-0.0124149646578933,-0.104396997365306,0.0122244874648869,0.0608515715420498,-0.0119125466162078,-0.0416075114003120,0.0114877299726374,0.0305139920646979,-0.0109610400592562,-0.0231767317862166,0.0103462462982944,0.0179123452201005,-0.00965872083006844,-0.0139392187955169,0.00891540962039221,0.0108448452085734,-0.00813404310905357,-0.00838972906356833,0.00733168780476571,0.00642411204097324,-0.00652616319850337,-0.00484707120349307,0.00573261054482327,0.00358652801599187,-0.00496649404133663,-0.00258731638098111,0.00423922975437062,0.00180683072389220,-0.00356275247645879,-0.00120602897645457,0.00294151183186525,0.000758811266874641,-0.00238769385139514,-0.000428686650077742,0.00189134412292593,0.000209215904988000,-0.00147731989791978,-4.63448770326659e-05,0.00109889241302976,-1.59234941706467e-05,-0.000838089356596170,0.000110049437588789,0.000530446929797254,-3.07332198143413e-05,-0.000469504941879963,0.000198295103017753,0.000330237662212286];
LP_base3 = [6.12780585326658e-05,1.63070602503443e-05,-6.72285628120998e-05,-0.000149541266354731,-0.000119810449222064,6.64308143667446e-05,0.000290030222885941,0.000316448110883590,8.04819102035179e-06,-0.000467862395439301,-0.000677624311760037,-0.000275486690927509,0.000583173950876424,0.00119767818321188,0.000841734299166718,-0.000490483308342665,-0.00181333942497140,-0.00180413276085515,-2.47189481995181e-05,0.00234918527705601,0.00318481888852270,0.00120393709872794,-0.00250892732918852,-0.00488029414334944,-0.00326687995154183,0.00187805480898982,0.00660955931045895,0.00634254676992654,4.60782051268906e-05,-0.00788122934317857,-0.0104061933867113,-0.00382125307915953,0.00796904879902335,0.0152381881139545,0.0100706124314083,-0.00585043479050356,-0.0204235286915247,-0.0196637672115522,-6.55765190337886e-05,0.0253991801027353,0.0345098061671373,0.0131558446933273,-0.0295443494098812,-0.0614042273667217,-0.0457651882292272,0.0322956747972992,0.150472416872618,0.257190402343197,0.300073519634591,0.257190402343197,0.150472416872618,0.0322956747972992,-0.0457651882292272,-0.0614042273667217,-0.0295443494098812,0.0131558446933273,0.0345098061671373,0.0253991801027353,-6.55765190337886e-05,-0.0196637672115522,-0.0204235286915247,-0.00585043479050356,0.0100706124314083,0.0152381881139545,0.00796904879902335,-0.00382125307915953,-0.0104061933867113,-0.00788122934317857,4.60782051268906e-05,0.00634254676992654,0.00660955931045895,0.00187805480898982,-0.00326687995154183,-0.00488029414334944,-0.00250892732918852,0.00120393709872794,0.00318481888852270,0.00234918527705601,-2.47189481995181e-05,-0.00180413276085515,-0.00181333942497140,-0.000490483308342665,0.000841734299166718,0.00119767818321188,0.000583173950876424,-0.000275486690927509,-0.000677624311760037,-0.000467862395439301,8.04819102035179e-06,0.000316448110883590,0.000290030222885941,6.64308143667446e-05,-0.000119810449222064,-0.000149541266354731,-6.72285628120998e-05,1.63070602503443e-05,6.12780585326658e-05];
% LP_base3 = [2.19038117042211e-05,-6.20841404157570e-05,-9.97611993011678e-05,-9.47729047964504e-05,-1.48464394967771e-06,0.000160291828657693,0.000291728363533567,0.000261011034362051,2.17091536635267e-06,-0.000394385200959140,-0.000681036159311772,-0.000582360677166308,-3.37159991487151e-06,0.000820098930855366,0.00137133401129982,0.00113913476811908,4.89276506094619e-06,-0.00153012594703318,-0.00250252363249786,-0.00203703106281935,-6.69283549991611e-06,0.00264595176452450,0.00425953936969232,0.00341739039927245,8.66960917944124e-06,-0.00433590433067623,-0.00690484726188547,-0.00548669845039959,-1.07358947179635e-05,0.00686378440593266,0.0108671149962161,0.00859759183657980,1.27283550754017e-05,-0.0107270333182705,-0.0169995751684161,-0.0134931481721584,-1.44957695423169e-05,0.0171112768435569,0.0274800858996100,0.0222176585053620,1.58989630767769e-05,-0.0299188511448005,-0.0503336563869308,-0.0434147896076745,-1.67997363136581e-05,0.0740368809674355,0.158230372677105,0.224763664690739,0.250017081102724,0.224763664690739,0.158230372677105,0.0740368809674355,-1.67997363136581e-05,-0.0434147896076745,-0.0503336563869308,-0.0299188511448005,1.58989630767769e-05,0.0222176585053620,0.0274800858996100,0.0171112768435569,-1.44957695423169e-05,-0.0134931481721584,-0.0169995751684161,-0.0107270333182705,1.27283550754017e-05,0.00859759183657980,0.0108671149962161,0.00686378440593266,-1.07358947179635e-05,-0.00548669845039959,-0.00690484726188547,-0.00433590433067623,8.66960917944124e-06,0.00341739039927245,0.00425953936969232,0.00264595176452450,-6.69283549991611e-06,-0.00203703106281935,-0.00250252363249786,-0.00153012594703318,4.89276506094619e-06,0.00113913476811908,0.00137133401129982,0.000820098930855366,-3.37159991487151e-06,-0.000582360677166308,-0.000681036159311772,-0.000394385200959140,2.17091536635267e-06,0.000261011034362051,0.000291728363533567,0.000160291828657693,-1.48464394967771e-06,-9.47729047964504e-05,-9.97611993011678e-05,-6.20841404157570e-05,2.19038117042211e-05];
% LP_base3 = [-6.68626527819765e-05,-7.22573663226062e-05,2.59724637286604e-05,0.000138030916179623,0.000168592864848942,-2.61428847141180e-05,-0.000287919073375483,-0.000349351394335339,4.79564571327912e-06,0.000517306474501090,0.000655462356510077,6.45013865431850e-05,-0.000844054338667943,-0.00113801062470633,-0.000221411104431505,0.00128171206232987,0.00185817967912811,0.000521883591014169,-0.00183663438288253,-0.00288786649249788,-0.00104173025330659,0.00250560643271562,0.00431199805603147,0.00188225440045845,-0.00327411454951619,-0.00623541268991753,-0.00318037345742469,0.00411578194481273,0.00879971145987546,0.00513103657753029,-0.00499325256444739,-0.0122228600123135,-0.00803783222769681,0.00586057705130150,0.0168949428962625,0.0124408360436902,-0.00666670934316737,-0.0236345561211203,-0.0194825341890430,0.00736016504313277,0.0345249696602514,0.0322258771032407,-0.00789401079610784,-0.0567555257151640,-0.0629175908662908,0.00823059610249638,0.141006537927771,0.271010448844496,0.324987733578038,0.271010448844496,0.141006537927771,0.00823059610249638,-0.0629175908662908,-0.0567555257151640,-0.00789401079610784,0.0322258771032407,0.0345249696602514,0.00736016504313277,-0.0194825341890430,-0.0236345561211203,-0.00666670934316737,0.0124408360436902,0.0168949428962625,0.00586057705130150,-0.00803783222769681,-0.0122228600123135,-0.00499325256444739,0.00513103657753029,0.00879971145987546,0.00411578194481273,-0.00318037345742469,-0.00623541268991753,-0.00327411454951619,0.00188225440045845,0.00431199805603147,0.00250560643271562,-0.00104173025330659,-0.00288786649249788,-0.00183663438288253,0.000521883591014169,0.00185817967912811,0.00128171206232987,-0.000221411104431505,-0.00113801062470633,-0.000844054338667943,6.45013865431850e-05,0.000655462356510077,0.000517306474501090,4.79564571327912e-06,-0.000349351394335339,-0.000287919073375483,-2.61428847141180e-05,0.000168592864848942,0.000138030916179623,2.59724637286604e-05,-7.22573663226062e-05,-6.68626527819765e-05];
% LP_base3 = [0.000351513380130003,0.000378039025224896,2.08270540755086e-05,-0.000218591742850680,-0.000633969891455429,-0.000514375783540978,-9.23796187788396e-05,0.000686594261580882,0.00110777490489252,0.000895974398309789,-0.000134854619102872,-0.00133217649169707,-0.00192255326300810,-0.00120790399717079,0.000578903268633281,0.00243610444381209,0.00296030193595920,0.00145444989792369,-0.00149571559331155,-0.00405051765836404,-0.00426274972205145,-0.00142794434416901,0.00306749612828441,0.00633988887436455,0.00575272152309883,0.000910061870622572,-0.00560658917282069,-0.00946487362500061,-0.00734990427667255,0.000462574966658700,0.00956359039240950,0.0137209327580219,0.00892949221851763,-0.00328700615682389,-0.0158110032850895,-0.0197626762219176,-0.0103562411959773,0.00885128743360873,0.0265255414323795,0.0295353226569302,0.0114939087210699,-0.0211411070292788,-0.0497150120648264,-0.0515807081072539,-0.0122280069412227,0.0650982228401988,0.158016085718602,0.233490222561327,0.262481645932529,0.233490222561327,0.158016085718602,0.0650982228401988,-0.0122280069412227,-0.0515807081072539,-0.0497150120648264,-0.0211411070292788,0.0114939087210699,0.0295353226569302,0.0265255414323795,0.00885128743360873,-0.0103562411959773,-0.0197626762219176,-0.0158110032850895,-0.00328700615682389,0.00892949221851763,0.0137209327580219,0.00956359039240950,0.000462574966658700,-0.00734990427667255,-0.00946487362500061,-0.00560658917282069,0.000910061870622572,0.00575272152309883,0.00633988887436455,0.00306749612828441,-0.00142794434416901,-0.00426274972205145,-0.00405051765836404,-0.00149571559331155,0.00145444989792369,0.00296030193595920,0.00243610444381209,0.000578903268633281,-0.00120790399717079,-0.00192255326300810,-0.00133217649169707,-0.000134854619102872,0.000895974398309789,0.00110777490489252,0.000686594261580882,-9.23796187788396e-05,-0.000514375783540978,-0.000633969891455429,-0.000218591742850680,2.08270540755086e-05,0.000378039025224896,0.000351513380130003];
CPM_recv = conv(CPM_recv, LP_base3);
CPM_recv = CPM_recv(48+1:end-48);
% figure
% plot((10*log10(abs(fft(CPM_recv)))))

% CPM_recv_one = CPM_recv(4:4:end);
% figure
% plot(angle(CPM_recv))
CPMHead = CPM_recv(4+1:15*4);
% CPMHeadDif = CPMHead(4+1:end) .*CPMHead(1:end-4);
CPMHeadAngle = 0;
for i = 1:7
    CPMHeadAngle = CPMHeadAngle + CPMHead((i-1)*4+1:4*i) .* conj(CPMHead((i-1+7)*4+1:4*(i+7)));
end

fre_est = angle(mean(CPMHeadAngle)) / 2 / pi / Ts / 7;

t = t(1:16:end);
CPM_recv = CPM_recv .* exp(-1i*2*pi*fd*t);

CPMHead = CPM_recv(4+1:15*4);

phi_est = angle(mean(CPMHead)) - pi/4;
% phi_est = mean(angle(CPMHead)) - pi/4;

CPM_recv_one = CPM_recv(15*4+1:end) .* exp(-1i*phi_est);
% CPM_recv_one = CPM_recv(15*4+1:end);
% plot(angle(CPM_recv))

% res = [out(1)];
% for i = 1:length(out)
%     if i > 1
%         res = [res, out(i-1)*out(i)];
%     end
% end
resCo = CPMdemodLikelyHead(CPM_recv_one);
% res = CPMdemod(CPM_recv_one);
% resD1 = CPMdemodD1(CPM_recv_one);
% resD2 = CPMdemodD2(CPM_recv_one);
% res(res == -1) = 0;
ACo = data - resCo;
ACo(end-1:end) = 0;
% AD1 = data - resD1;
% AD1(end-1:end) = 0;
% AD2 = data - resD2;
% AD2(end-1:end) = 0;

erro_cnt_Co = erro_cnt_Co + sum(abs(ACo));
erro_cnt_D1 = erro_cnt_D1 + sum(abs(AD1));
erro_cnt_D2 = erro_cnt_D2 + sum(abs(AD2));
end
Co_rec = [Co_rec, erro_cnt_Co / Ne / (48-15)];
D1_rec = [D1_rec, erro_cnt_D1 / Ne / (48-15)];
D2_rec = [D2_rec, erro_cnt_D2 / Ne / (48-15)];
end

ebn0 = 0:10;
f = figure;
f.PaperUnits = 'centimeters';
f.PaperSize = [16, 12];
f.Units = 'centimeters';
f.Position = [0, 0, 16, 12];
semilogy(ebn0, Co_rec, '-s', 'LineWidth', 2);
hold on;
semilogy(ebn0, D1_rec, '-s', 'LineWidth', 2);
hold on;
semilogy(ebn0, D2_rec, '-s', 'LineWidth', 2);
hold on;
grid on;

legend('相干viterbi','D1','D2')
