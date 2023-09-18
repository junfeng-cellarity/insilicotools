package com.insilico.application.insilicotools.gui.filter.eganEgg;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterState;


import java.awt.*;
import java.awt.geom.Path2D;
import java.util.ArrayList;

public class EganEggFilter extends Filter {
    private final Shape abs95Shape;
    private final Shape abs99Shape;
    private final Shape bbb90Shape;
    private final Shape bbb99Shape;

    private final double[][] ABS95Pct = {{129.996535, 0.486020186}, {128.6612167, -0.292001894}, {127.3374055, -0.589993822}, {126.015528, -0.807320692}, {124.6946677, -0.982213833}, {123.3744679, -1.129555401}, {122.0547447, -1.257012515}, {120.7353885, -1.369159108}, {119.4163277, -1.468987181}, {118.0975121, -1.558582677}, {116.7789053, -1.639469554}, {115.4604796, -1.712801878}, {114.1422134, -1.779478766}, {112.8240896, -1.840217014}, {111.5060941, -1.895599026}, {110.1882156, -1.94610562}, {108.8704443, -1.992139153}, {107.552772, -2.034040261}, {106.2351918, -2.07210024}, {104.9176977, -2.106570383}, {103.6002844, -2.137669142}, {102.2829473, -2.165587705}, {100.9656824, -2.190494392}, {99.64848615, -2.212538174}, {98.33135535, -2.231851502}, {97.01428717, -2.248552622}, {95.69727907, -2.262747474}, {94.38032876, -2.274531264}, {93.0634342, -2.283989779}, {91.74659352, -2.291200494}, {90.42980504, -2.296233501}, {89.11306725, -2.299152303}, {87.79637875, -2.300014492}, {86.47973831, -2.29887233}, {85.16314477, -2.295773249}, {83.84659711, -2.290760278}, {82.5300944, -2.283872424}, {81.21363578, -2.275144991}, {79.8972205, -2.264609863}, {78.58084786, -2.25229575}, {77.26451725, -2.2382284}, {75.94822812, -2.222430785}, {74.63197998, -2.204923258}, {73.3157724, -2.185723695}, {71.999605, -2.164847612}, {70.68347748, -2.142308265}, {69.36738956, -2.118116735}, {68.05134104, -2.092281995}, {66.73533174, -2.064810972}, {65.41936154, -2.03570858}, {64.10343038, -2.004977758}, {62.78753824, -1.972619482}, {61.47168513, -1.938632777}, {60.15587112, -1.903014703}, {58.84009634, -1.865760345}, {57.52436096, -1.826862779}, {56.20866517, -1.786313028}, {54.89300926, -1.744100014}, {53.57739354, -1.700210479}, {52.26181838, -1.654628909}, {50.9462842, -1.607337427}, {49.63079151, -1.558315679}, {48.31534084, -1.507540695}, {46.99993282, -1.454986725}, {45.68456813, -1.400625061}, {44.36924754, -1.344423818}, {43.05397189, -1.286347691}, {41.73874212, -1.226357676}, {40.42355926, -1.16441074}, {39.10842445, -1.100459453}, {37.79333894, -1.034451554}, {36.47830411, -0.96632945}, {35.16332149, -0.896029638}, {33.84839275, -0.823482025}, {32.53351974, -0.748609138}, {31.21870454, -0.671325189}, {29.90394942, -0.591534971}, {28.58925691, -0.509132546}, {27.27462986, -0.423999667}, {25.96007143, -0.336003882}, {24.64558521, -0.244996222}, {23.33117522, -0.150808366}, {22.01684604, -0.053249125}, {20.70260293, 0.047899951}, {19.38845189, 0.152890157}, {18.07439991, 0.262012787}, {16.76045515, 0.375608478}, {15.44662727, 0.494079587}, {14.13292779, 0.617906933}, {12.81937066, 0.747672919}, {11.50597308, 0.88409434}, {10.19275659, 1.028070314}, {8.879748866, 1.180754907}, {7.566986429, 1.343672078}, {6.254519288, 1.518907768}, {4.942419164, 1.709453979}, {3.630795704, 1.919884645}, {2.319832699, 2.157866883}, {1.009886901, 2.43828285}, {-0.298125225, 2.799363876}, {-1.594630245, 3.640475052}, {-1.594630245, 3.640475052}, {-0.25931196, 4.418497132}, {1.064499219, 4.716489061}, {2.386376726, 4.93381593}, {3.707237026, 5.108709071}, {5.027436871, 5.256050639}, {6.347160053, 5.383507753}, {7.666516216, 5.495654346}, {8.985577085, 5.595482419}, {10.30439266, 5.685077915}, {11.62299948, 5.765964792}, {12.9414252, 5.839297116}, {14.25969138, 5.905974004}, {15.5778152, 5.966712252}, {16.89581063, 6.022094265}, {18.21368918, 6.072600858}, {19.53146051, 6.118634391}, {20.84913277, 6.160535499}, {22.16671296, 6.198595478}, {23.48420709, 6.233065621}, {24.8016204, 6.26416438}, {26.11895748, 6.292082943}, {27.43622236, 6.31698963}, {28.75341862, 6.339033412}, {30.07054942, 6.35834674}, {31.3876176, 6.37504786}, {32.7046257, 6.389242712}, {34.021576, 6.401026502}, {35.33847056, 6.410485017}, {36.65531124, 6.417695732}, {37.97209972, 6.422728739}, {39.28883751, 6.425647541}, {40.60552601, 6.42650973}, {41.92216646, 6.425367568}, {43.23875999, 6.422268487}, {44.55530765, 6.417255516}, {45.87181036, 6.410367662}, {47.18826898, 6.401640229}, {48.50468426, 6.391105101}, {49.8210569, 6.378790988}, {51.13738751, 6.364723638}, {52.45367664, 6.348926023}, {53.76992478, 6.331418496}, {55.08613237, 6.312218933}, {56.40229976, 6.29134285}, {57.71842728, 6.268803503}, {59.0345152, 6.244611973}, {60.35056372, 6.218777234}, {61.66657302, 6.19130621}, {62.98254322, 6.162203818}, {64.29847438, 6.131472996}, {65.61436652, 6.099114721}, {66.93021964, 6.065128015}, {68.24603364, 6.029509942}, {69.56180842, 5.992255584}, {70.87754381, 5.953358017}, {72.19323959, 5.912808266}, {73.5088955, 5.870595252}, {74.82451122, 5.826705717}, {76.14008639, 5.781124147}, {77.45562056, 5.733832665}, {78.77111326, 5.684810917}, {80.08656392, 5.634035933}, {81.40197195, 5.581481963}, {82.71733663, 5.527120299}, {84.03265723, 5.470919056}, {85.34793287, 5.412842929}, {86.66316264, 5.352852914}, {87.9783455, 5.290905978}, {89.29348031, 5.226954691}, {90.60856582, 5.160946792}, {91.92360065, 5.092824688}, {93.23858327, 5.022524876}, {94.55351202, 4.949977263}, {95.86838502, 4.875104376}, {97.18320022, 4.797820427}, {98.49795535, 4.718030209}, {99.81264785, 4.635627784}, {101.1272749, 4.550494905}, {102.4418333, 4.46249912}, {103.7563196, 4.37149146}, {105.0707295, 4.277303604}, {106.3850587, 4.179744363}, {107.6993018, 4.078595287}, {109.0134529, 3.973605082}, {110.3275049, 3.864482452}, {111.6414496, 3.75088676}, {112.9552775, 3.632415651}, {114.268977, 3.508588306}, {115.5825341, 3.37882232}, {116.8959317, 3.242400898}, {118.2091482, 3.098424924}, {119.5221559, 2.945740331}, {120.8349183, 2.78282316}, {122.1473855, 2.60758747}, {123.4594856, 2.417041259}, {124.7711091, 2.206610593}, {126.0820721, 1.968628355}, {127.3920179, 1.688212388}, {128.70003, 1.327131363}, {129.996535, 0.486020186}};
    private final double[][] ABS99Pct = {{146.0961678, 0.100085728}, {144.4341082, -0.868311894}, {142.7863714, -1.239219973}, {141.1410415, -1.509724923}, {139.4969776, -1.727412964}, {137.8537358, -1.91080779}, {136.2110873, -2.069452601}, {134.5688957, -2.209040535}, {132.9270716, -2.333295709}, {131.2855528, -2.444814478}, {129.6442938, -2.545493702}, {128.0032603, -2.636769837}, {126.3624253, -2.719762007}, {124.7217675, -2.795362399}, {123.0812695, -2.864295931}, {121.4409171, -2.927161068}, {119.800698, -2.984458624}, {118.1606023, -3.036612585}, {116.5206212, -3.083985525}, {114.8807472, -3.126890226}, {113.2409738, -3.165598594}, {111.6012952, -3.200348597}, {109.9617066, -3.231349744}, {108.3222033, -3.258787456}, {106.6827816, -3.282826595}, {105.0434377, -3.30361434}, {103.4041687, -3.321282555}, {101.7649716, -3.33594974}, {100.1258439, -3.347722676}, {98.4867832, -3.356697793}, {96.84778751, -3.362962334}, {95.20885491, -3.366595343}, {93.56998368, -3.367668503}, {91.93117224, -3.366246864}, {90.2924192, -3.362389462}, {88.65372326, -3.356149859}, {87.01508326, -3.347576604}, {85.37649816, -3.336713641}, {83.73796698, -3.323600654}, {82.09948889, -3.308273381}, {80.46106311, -3.290763867}, {78.82268896, -3.271100706}, {77.18436582, -3.249309233}, {75.54609317, -3.225411697}, {73.90787055, -3.19942741}, {72.26969755, -3.171372872}, {70.63157385, -3.141261875}, {68.99349917, -3.109105589}, {67.35547332, -3.074912634}, {65.71749615, -3.038689128}, {64.07956756, -3.000438729}, {62.44168753, -2.960162653}, {60.80385609, -2.917859683}, {59.16607333, -2.873526162}, {57.52833938, -2.827155972}, {55.89065447, -2.778740494}, {54.25301885, -2.728268556}, {52.61543286, -2.675726367}, {50.9778969, -2.621097428}, {49.34041141, -2.564362426}, {47.70297695, -2.505499111}, {46.06559412, -2.444482149}, {44.42826359, -2.381282947}, {42.79098615, -2.315869458}, {41.15376264, -2.248205946}, {39.51659402, -2.178252725}, {37.87948135, -2.105965852}, {36.24242577, -2.031296778}, {34.60542859, -1.954191941}, {32.96849121, -1.874592306}, {31.3316152, -1.792432822}, {29.69480226, -1.707641805}, {28.05805431, -1.620140212}, {26.42137342, -1.5298408}, {24.78476191, -1.43664714}, {23.14822234, -1.34045245}, {21.51175756, -1.241138228}, {19.87537071, -1.138572614}, {18.23906533, -1.032608426}, {16.60284538, -0.923080803}, {14.96671529, -0.809804325}, {13.33068008, -0.692569482}, {11.69474548, -0.571138305}, {10.05891798, -0.44523889}, {8.42320509, -0.314558453}, {6.787615502, -0.178734421}, {5.152159378, -0.037342808}, {3.516848723, 0.110117199}, {1.881697883, 0.264244067}, {0.246724237, 0.425762712}, {-1.388050829, 0.595565323}, {-3.022600488, 0.774771023}, {-4.656890305, 0.964816268}, {-6.29087481, 1.167597917}, {-7.924491764, 1.385712326}, {-9.557651893, 1.622883613}, {-11.19021872, 1.884804915}, {-12.82196349, 2.181019431}, {-14.45244215, 2.530050857}, {-16.08051398, 2.979485412}, {-17.69426302, 4.02640951}, {-17.69426302, 4.02640951}, {-16.03220343, 4.994807132}, {-14.38446665, 5.365715212}, {-12.73913669, 5.636220161}, {-11.09507284, 5.853908202}, {-9.451831059, 6.037303028}, {-7.809182573, 6.19594784}, {-6.166990911, 6.335535773}, {-4.5251668, 6.459790947}, {-2.883648001, 6.571309716}, {-1.242389044, 6.67198894}, {0.398644506, 6.763265075}, {2.039479475, 6.846257245}, {3.680137252, 6.921857637}, {5.320635213, 6.990791169}, {6.960987704, 7.053656306}, {8.601206732, 7.110953862}, {10.24130246, 7.163107823}, {11.88128358, 7.210480763}, {13.52115759, 7.253385464}, {15.160931, 7.292093832}, {16.80060953, 7.326843835}, {18.44019819, 7.357844982}, {20.07970142, 7.385282694}, {21.71912319, 7.409321833}, {23.35846702, 7.430109579}, {24.99773607, 7.447777793}, {26.63693318, 7.462444978}, {28.27606091, 7.474217914}, {29.91512157, 7.483193031}, {31.55411725, 7.489457572}, {33.19304985, 7.493090581}, {34.83192109, 7.494163741}, {36.47073252, 7.492742102}, {38.10948556, 7.4888847}, {39.7481815, 7.482645097}, {41.3868215, 7.474071843}, {43.02540661, 7.463208879}, {44.66393778, 7.450095893}, {46.30241587, 7.434768619}, {47.94084165, 7.417259105}, {49.5792158, 7.397595944}, {51.21753894, 7.375804471}, {52.85581159, 7.351906935}, {54.49403421, 7.325922648}, {56.13220721, 7.29786811}, {57.77033092, 7.267757113}, {59.40840559, 7.235600827}, {61.04643144, 7.201407872}, {62.68440861, 7.165184366}, {64.3223372, 7.126933967}, {65.96021723, 7.086657891}, {67.59804867, 7.044354921}, {69.23583144, 7.0000214}, {70.87356538, 6.95365121}, {72.51125029, 6.905235732}, {74.14888591, 6.854763794}, {75.7864719, 6.802221605}, {77.42400787, 6.747592666}, {79.06149335, 6.690857664}, {80.69892781, 6.631994349}, {82.33631065, 6.570977387}, {83.97364117, 6.507778185}, {85.61091861, 6.442364696}, {87.24814212, 6.374701184}, {88.88531074, 6.304747963}, {90.52242342, 6.23246109}, {92.15947899, 6.157792016}, {93.79647617, 6.080687179}, {95.43341355, 6.001087544}, {97.07028957, 5.91892806}, {98.7071025, 5.834137043}, {100.3438505, 5.74663545}, {101.9805313, 5.656336038}, {103.6171428, 5.563142378}, {105.2536824, 5.466947688}, {106.8901472, 5.367633466}, {108.5265341, 5.265067852}, {110.1628394, 5.159103664}, {111.7990594, 5.049576042}, {113.4351895, 4.936299563}, {115.0712247, 4.81906472}, {116.7071593, 4.697633543}, {118.3429868, 4.571734128}, {119.9786997, 4.441053691}, {121.6142893, 4.305229659}, {123.2497454, 4.163838046}, {124.885056, 4.016378039}, {126.5202069, 3.862251171}, {128.1551805, 3.700732526}, {129.7899556, 3.530929915}, {131.4245052, 3.351724215}, {133.0587951, 3.16167897}, {134.6927796, 2.958897321}, {136.3263965, 2.740782912}, {137.9595567, 2.503611625}, {139.5921235, 2.241690323}, {141.2238683, 1.945475807}, {142.8543469, 1.596444381}, {144.4824187, 1.147009827}, {146.0961678, 0.100085728}};
    private final double[][] BBB90Pct = {{-5.236517058, 5.076807435}, {-4.243907377, 5.624621148}, {-3.265457974, 5.824249954}, {-2.289388088, 5.965369193}, {-1.314569943, 6.075709566}, {-0.340564535, 6.16606569}, {0.632854308, 6.241998831}, {1.60582151, 6.306826644}, {2.578425331, 6.362719346}, {3.550727304, 6.411189953}, {4.522772383, 6.453343854}, {5.494594613, 6.490018138}, {6.466220517, 6.521864971}, {7.437671237, 6.549404272}, {8.408963956, 6.573058483}, {9.380112856, 6.59317636}, {10.35112981, 6.610049748}, {11.32202486, 6.623925724}, {12.2928066, 6.635015577}, {13.26348244, 6.64350158}, {14.23405884, 6.649542185}, {15.20454142, 6.653276066}, {16.17493515, 6.654825315}, {17.14524444, 6.654297985}, {18.11547317, 6.651790149}, {19.08562486, 6.647387575}, {20.05570261, 6.641167105}, {21.02570923, 6.633197796}, {21.99564727, 6.623541873}, {22.96551899, 6.612255531}, {23.93532648, 6.599389612}, {24.9050716, 6.584990175}, {25.87475605, 6.569098996}, {26.84438138, 6.551753982}, {27.81394898, 6.532989536}, {28.78346013, 6.51283687}, {29.75291596, 6.491324277}, {30.72231753, 6.468477365}, {31.69166578, 6.44431926}, {32.66096155, 6.418870786}, {33.6302056, 6.39215062}, {34.59939861, 6.364175424}, {35.56854118, 6.334959961}, {36.53763384, 6.304517198}, {37.50667704, 6.272858388}, {38.47567118, 6.239993146}, {39.44461658, 6.205929511}, {40.41351351, 6.17067399}, {41.38236217, 6.134231608}, {42.3511627, 6.09660593}, {43.3199152, 6.057799086}, {44.28861969, 6.017811787}, {45.25727615, 5.976643322}, {46.22588448, 5.934291562}, {47.19444454, 5.890752939}, {48.16295614, 5.846022432}, {49.13141899, 5.800093531}, {50.09983278, 5.752958198}, {51.06819711, 5.704606818}, {52.03651153, 5.655028138}, {53.00477552, 5.604209192}, {53.97298846, 5.552135215}, {54.94114968, 5.498789546}, {55.90925842, 5.444153508}, {56.87731384, 5.388206278}, {57.845315, 5.330924728}, {58.81326085, 5.272283251}, {59.78115024, 5.212253554}, {60.7489819, 5.150804425}, {61.71675444, 5.087901461}, {62.68446631, 5.023506755}, {63.65211581, 4.957578532}, {64.61970108, 4.890070731}, {65.58722004, 4.82093251}, {66.5546704, 4.750107677}, {67.52204964, 4.677534004}, {68.48935495, 4.603142435}, {69.45658321, 4.526856129}, {70.42373091, 4.448589316}, {71.39079417, 4.368245924}, {72.35776858, 4.285717899}, {73.32464918, 4.200883152}, {74.29143032, 4.113603006}, {75.25810557, 4.02371901}, {76.22466752, 3.931048891}, {77.19110756, 3.83538136}, {78.15741565, 3.736469341}, {79.12357992, 3.634020987}, {80.08958619, 3.527687543}, {81.05541728, 3.417046568}, {82.02105204, 3.301578141}, {82.98646395, 3.180630097}, {83.95161897, 3.053365347}, {84.91647214, 2.918678502}, {85.88096193, 2.775056546}, {86.84500008, 2.620329262}, {87.80845166, 2.451178995}, {88.77109051, 2.26204448}, {89.73247762, 2.042131098}, {90.6914852, 1.763708149}, {91.63633251, 1.137100293}, {90.64372283, 0.589286579}, {89.66527343, 0.389657773}, {88.68920354, 0.248538534}, {87.7143854, 0.138198162}, {86.74037999, 0.047842038}, {85.76696115, -0.028091103}, {84.79399395, -0.092918916}, {83.82139013, -0.148811618}, {82.84908815, -0.197282225}, {81.87704307, -0.239436126}, {80.90522084, -0.27611041}, {79.93359494, -0.307957243}, {78.96214422, -0.335496544}, {77.9908515, -0.359150755}, {77.0197026, -0.379268632}, {76.04868565, -0.39614202}, {75.0777906, -0.410017996}, {74.10700886, -0.421107849}, {73.13633301, -0.429593852}, {72.16575662, -0.435634457}, {71.19527404, -0.439368338}, {70.2248803, -0.440917587}, {69.25457102, -0.440390257}, {68.28434228, -0.437882421}, {67.3141906, -0.433479847}, {66.34411285, -0.427259377}, {65.37410622, -0.419290068}, {64.40416819, -0.409634145}, {63.43429646, -0.398347804}, {62.46448898, -0.385481884}, {61.49474385, -0.371082447}, {60.5250594, -0.355191268}, {59.55543407, -0.337846254}, {58.58586647, -0.319081808}, {57.61635533, -0.298929143}, {56.64689949, -0.27741655}, {55.67749792, -0.254569637}, {54.70814967, -0.230411532}, {53.73885391, -0.204963058}, {52.76960985, -0.178242892}, {51.80041684, -0.150267696}, {50.83127427, -0.121052233}, {49.86218162, -0.09060947}, {48.89313841, -0.05895066}, {47.92414428, -0.026085418}, {46.95519888, 0.007978217}, {45.98630195, 0.043233737}, {45.01745329, 0.07967612}, {44.04865275, 0.117301798}, {43.07990025, 0.156108641}, {42.11119576, 0.196095941}, {41.14253931, 0.237264405}, {40.17393097, 0.279616166}, {39.20537091, 0.323154789}, {38.23685932, 0.367885296}, {37.26839647, 0.413814197}, {36.29998268, 0.46094953}, {35.33161834, 0.509300909}, {34.36330392, 0.558879589}, {33.39503994, 0.609698536}, {32.426827, 0.661772513}, {31.45866578, 0.715118182}, {30.49055703, 0.76975422}, {29.52250161, 0.82570145}, {28.55450046, 0.882983}, {27.58655461, 0.941624477}, {26.61866522, 1.001654174}, {25.65083356, 1.063103303}, {24.68306102, 1.126006267}, {23.71534915, 1.190400973}, {22.74769964, 1.256329196}, {21.78011438, 1.323836997}, {20.81259542, 1.392975218}, {19.84514505, 1.463800051}, {18.87776581, 1.536373724}, {17.9104605, 1.610765293}, {16.94323225, 1.687051599}, {15.97608454, 1.765318412}, {15.00902129, 1.845661804}, {14.04204688, 1.928189828}, {13.07516628, 2.013024576}, {12.10838513, 2.100304721}, {11.14170988, 2.190188717}, {10.17514794, 2.282858837}, {9.208707899, 2.378526367}, {8.242399807, 2.477438387}, {7.276235534, 2.579886741}, {6.310229264, 2.686220185}, {5.344398176, 2.79686116}, {4.378763414, 2.912329587}, {3.413351503, 3.033277631}, {2.448196484, 3.160542381}, {1.483343314, 3.295229226}, {0.518853524, 3.438851182}, {-0.445184625, 3.593578465}, {-1.408636208, 3.762728732}, {-2.371275054, 3.951863248}, {-3.33266216, 4.17177663}, {-4.291669749, 4.450199578}, {-5.236517058, 5.076807435}};
    private final double[][] BBB99Pct = {{-25.39089431, 5.896462729}, {-23.98526014, 6.672221513}, {-22.59967832, 6.954915786}, {-21.21746614, 7.154754685}, {-19.83702655, 7.311007643}, {-18.45773788, 7.438960915}, {-17.07927984, 7.546489806}, {-15.70146137, 7.638292446}, {-14.32415748, 7.717442079}, {-12.94728104, 7.786081287}, {-11.57076838, 7.84577541}, {-10.19457131, 7.897709849}, {-8.818652248, 7.942808136}, {-7.442981265, 7.981806529}, {-6.067534029, 8.015303249}, {-4.692290455, 8.043792165}, {-3.317233734, 8.067686564}, {-1.942349638, 8.087336328}, {-0.567625998, 8.103040665}, {0.806947684, 8.11505769}, {2.181380532, 8.123611788}, {3.555680534, 8.128899336}, {4.92985472, 8.131093227}, {6.303909313, 8.130346475}, {7.677849847, 8.126795129}, {9.051681261, 8.120560646}, {10.42540798, 8.111751841}, {11.79903398, 8.100466505}, {13.17256285, 8.086792757}, {14.54599782, 8.070810173}, {15.91934182, 8.052590749}, {17.2925975, 8.032199713}, {18.66576727, 8.009696221}, {20.03885332, 7.985133953}, {21.41185762, 7.958561628}, {22.78478196, 7.930023446}, {24.15762799, 7.899559472}, {25.53039717, 7.867205968}, {26.90309085, 7.832995685}, {28.2757102, 7.796958111}, {29.64825632, 7.759119695}, {31.02073016, 7.719504031}, {32.39313257, 7.678132025}, {33.7654643, 7.63502204}, {35.137726, 7.590190012}, {36.50991821, 7.543649556}, {37.88204141, 7.495412055}, {39.25409597, 7.445486727}, {40.62608217, 7.393880684}, {41.99800023, 7.340598976}, {43.36985026, 7.28564462}, {44.74163231, 7.229018621}, {46.11334634, 7.170719975}, {47.48499221, 7.110745663}, {48.85656974, 7.049090637}, {50.22807862, 6.985747783}, {51.59951849, 6.920707884}, {52.97088887, 6.853959558}, {54.34218922, 6.785489189}, {55.71341889, 6.715280839}, {57.08457714, 6.643316149}, {58.4556631, 6.569574211}, {59.82667583, 6.49403143}, {61.19761424, 6.416661358}, {62.56847713, 6.337434508}, {63.93926319, 6.256318127}, {65.30997092, 6.173275955}, {66.68059871, 6.088267925}, {68.05114474, 6.001249838}, {69.42160705, 5.912172975}, {70.79198345, 5.820983657}, {72.16227153, 5.727622726}, {73.53246865, 5.632024956}, {74.90257186, 5.534118349}, {76.27257794, 5.43382333}, {77.6424833, 5.331051781}, {79.01228397, 5.22570591}, {80.38197552, 5.117676901}, {81.75155301, 5.006843298}, {83.1210109, 4.893069053}, {84.49034298, 4.77620115}, {85.85954222, 4.656066698}, {87.22860062, 4.532469318}, {88.59750906, 4.405184626}, {89.96625704, 4.273954507}, {91.3348324, 4.138479754}, {92.70322091, 3.998410483}, {94.07140576, 3.853333409}, {95.43936686, 3.702754661}, {96.80707988, 3.546076019}, {98.17451488, 3.382561226}, {99.54163431, 3.211286748}, {100.90839, 3.031067186}, {102.2747181, 2.840337198}, {103.6405318, 2.636954203}, {105.0057058, 2.417844958}, {106.3700492, 2.178311331}, {107.7332417, 1.910478019}, {109.0946616, 1.599058765}, {110.4527119, 1.204784137}, {111.7907098, 0.317444999}, {110.3850756, -0.458313785}, {108.9994938, -0.741008058}, {107.6172816, -0.940846957}, {106.236842, -1.097099915}, {104.8575533, -1.225053187}, {103.4790953, -1.332582078}, {102.1012768, -1.424384718}, {100.7239729, -1.503534351}, {99.34709649, -1.572173559}, {97.97058384, -1.631867682}, {96.59438676, -1.683802121}, {95.2184677, -1.728900408}, {93.84279672, -1.767898802}, {92.46734949, -1.801395521}, {91.09210591, -1.829884438}, {89.71704919, -1.853778836}, {88.34216509, -1.8734286}, {86.96744145, -1.889132937}, {85.59286777, -1.901149963}, {84.21843492, -1.90970406}, {82.84413492, -1.914991608}, {81.46996074, -1.917185499}, {80.09590614, -1.916438748}, {78.72196561, -1.912887402}, {77.34813419, -1.906652918}, {75.97440748, -1.897844113}, {74.60078147, -1.886558777}, {73.22725261, -1.872885029}, {71.85381764, -1.856902445}, {70.48047364, -1.838683021}, {69.10721795, -1.818291985}, {67.73404818, -1.795788493}, {66.36096214, -1.771226225}, {64.98795784, -1.7446539}, {63.61503349, -1.716115719}, {62.24218746, -1.685651744}, {60.86941828, -1.65329824}, {59.49672461, -1.619087957}, {58.12410525, -1.583050384}, {56.75155914, -1.545211967}, {55.3790853, -1.505596303}, {54.00668288, -1.464224298}, {52.63435116, -1.421114312}, {51.26208946, -1.376282284}, {49.88989724, -1.329741829}, {48.51777404, -1.281504327}, {47.14571949, -1.231578999}, {45.77373328, -1.179972956}, {44.40181522, -1.126691248}, {43.02996519, -1.071736892}, {41.65818314, -1.015110893}, {40.28646912, -0.956812247}, {38.91482324, -0.896837935}, {37.54324572, -0.835182909}, {36.17173683, -0.771840056}, {34.80029697, -0.706800156}, {33.42892658, -0.64005183}, {32.05762623, -0.571581461}, {30.68639656, -0.501373112}, {29.31523832, -0.429408421}, {27.94415236, -0.355666483}, {26.57313963, -0.280123702}, {25.20220122, -0.20275363}, {23.83133832, -0.12352678}, {22.46055227, -0.0424104}, {21.08984453, 0.040631773}, {19.71921675, 0.125639803}, {18.34867071, 0.21265789}, {16.9782084, 0.301734752}, {15.607832, 0.392924071}, {14.23754392, 0.486285001}, {12.86734681, 0.581882772}, {11.4972436, 0.679789379}, {10.12723752, 0.780084398}, {8.757332153, 0.882855947}, {7.387531485, 0.988201818}, {6.017839937, 1.096230827}, {4.648262449, 1.20706443}, {3.278804553, 1.320838675}, {1.909472473, 1.437706578}, {0.54027324, 1.55784103}, {-0.82878516, 1.68143841}, {-2.197693602, 1.808723102}, {-3.566441587, 1.939953221}, {-4.935016948, 2.075427973}, {-6.303405455, 2.215497245}, {-7.671590301, 2.360574319}, {-9.0395514, 2.511153067}, {-10.40726442, 2.667831709}, {-11.77469943, 2.831346502}, {-13.14181886, 3.00262098}, {-14.5085745, 3.182840542}, {-15.87490269, 3.37357053}, {-17.2407163, 3.576953525}, {-18.60589034, 3.79606277}, {-19.97023375, 4.035596396}, {-21.33342624, 4.303429709}, {-22.69484614, 4.614848963}, {-24.05289641, 5.00912359}, {-25.39089431, 5.896462729}};


    public EganEggFilter() {
        abs95Shape = createPolygon(ABS95Pct);
        abs99Shape = createPolygon(ABS99Pct);
        bbb90Shape = createPolygon(BBB90Pct);
        bbb99Shape = createPolygon(BBB99Pct);
    }

    private Shape createPolygon(double[][] points) {
        Path2D polygon = new Path2D.Double();
        for (int i = 0; i < points.length; i++) {
        	if (i==0) 
        		polygon.moveTo(points[i][0], points[i][1]);
        	else 
        		polygon.lineTo(points[i][0], points[i][1]);

        }
        return polygon;
    }

    public FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] molecules, final FilterState state) throws Exception {
        final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
        final ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();

        final EganEggFilterState filterState = (EganEggFilterState) state;

        if (molecules.length > 0) {

            for(PropertyMolecule m:molecules) {
                double psa = m.getPSA();
                double clogp = m.getCLogP();
                if (filterState.getAbsState() == EganEggFilterState.IGNORE || (filterState.getAbsState() == EganEggFilterState.ABS_95 && abs95Shape.contains(psa, clogp)) || (filterState.getAbsState() == EganEggFilterState.ABS_99 && abs99Shape.contains(psa, clogp))) {
                    if (filterState.getBbbState() == EganEggFilterState.IGNORE || (filterState.getBbbState() == EganEggFilterState.BBB_90 && bbb90Shape.contains(psa, clogp)) || (filterState.getBbbState() == EganEggFilterState.BBB_99 && bbb99Shape.contains(psa, clogp))) {
                        pass.add(m);
                    } else {
                        fail.add(m);
                    }
                } else {
                    fail.add(m);
                }
            }

        }
        return new FilterResult(pass.toArray(new PropertyMolecule[pass.size()]), fail.toArray(new PropertyMolecule[fail.size()]));
    }
}