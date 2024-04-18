from pathlib import Path

GOLDEN_DATA_PATH = Path(Path(__file__).parent, "golden_data")


PAM_EXPECTED_SCORES = [
    0.001239810431396115,
    0.006321719197442959,
    0.022366417671350794,
    0.02356005692376431,
    0.023591089801209657,
    0.0237015155252253,
    0.02494871077780501,
    0.02765531021502237,
    0.05723352359279142,
    -0.026031899697247438,
    -0.045399225810167206,
    -0.048802117011236074,
    -0.049071670232245686,
    -0.050609092352903964,
    -0.052005958261434004,
    -0.05212809666055347,
    -0.05219436481287494,
    -0.05229549647527615,
    -0.05562289106810496,
    -0.0680586795708448,
    -0.080486827060348,
    -0.08708614518433638,
    -0.09347270442733403,
    -0.10733976984019238,
    -0.09525599973860271,
    -0.09344614008515476,
    -0.09203987891656933,
    -0.07762954905267365,
    -0.08947079582415686,
    -0.09173525480186309,
    -0.09743570341595072,
    -0.10060734372517446,
    -0.13928893297991315,
    -0.17460375352669397,
    -0.17475261110396872,
    -0.17031599328406574,
    -0.1661587701497279,
    -0.1664421078509678,
    -0.16478663784890032,
    -0.16609846719661853,
    -0.16691114907024981,
    -0.1673165715851486,
    -0.17164095778004496,
    -0.1732524510519868,
    -0.17352716673977173,
    -0.17394908243845958,
    -0.17407532356221084,
    -0.1697432597263014,
    -0.17090737500355396,
    -0.17148612056561355,
    -0.17243034254792916,
    -0.17250242028599208,
    -0.17349744530607505,
    -0.17156916314048715,
    -0.17087689628648867,
    -0.17039645590656666,
    -0.17124201869844918,
    -0.17122651284437534,
    -0.17238041512299634,
    -0.17269067870892135,
    -0.1738149625061334,
    -0.17480998867178196,
    -0.17479150234600335,
    -0.17488886861531663,
    -0.17514401312700567,
    -0.17517583629278174,
    -0.1753688433367886,
    -0.17640234267324945,
    -0.1772456623081938,
    -0.17716139852885898,
    -0.17630326766858379,
    -0.17673989971707982,
    -0.17734885123070754,
    -0.17721649370316805,
    -0.1777042024199318,
    -0.17826981243725237,
    -0.17854430105516525,
    -0.17856624282100253,
    -0.17869271946472667,
    -0.17875479287980137,
    -0.17878826358419384,
    -0.17970716028636985,
    -0.1806352773414821,
    -0.18023235943323299,
    -0.1800376011216546,
    -0.18075378827525296,
    -0.18163713539039433,
    -0.18234358338205933,
    -0.1814940376576018,
    -0.18161785274664427,
    -0.1816573858126205,
    -0.1821612189289432,
    -0.1815902086414353,
    -0.18115171590553097,
    -0.18117653831295877,
    -0.1811891114798373,
    -0.18131709218972805,
    -0.18067410292225058,
    -0.18104446672682117,
    -0.18106004855540903,
    -0.18107784463968024,
    -0.18137076403843633,
    -0.17984872596031024,
    -0.17992682094263204,
    -0.17995498721146755,
    -0.17997386557447312,
    -0.18004579052015068,
    -0.18006692743072406,
    -0.18032752189520404,
    -0.18034793878068686,
    -0.180549332168284,
    -0.18057405694667325,
    -0.1807566033056329,
    -0.18093253195751252,
    -0.18096579039053484,
    -0.18098671009042228,
    -0.18157643652459057,
    -0.18181564752199705,
    -0.18228752384764443,
    -0.18178799062428616,
    -0.18179150375381325,
    -0.1818562866740624,
    -0.1818842938922417,
    -0.18190309498137486,
    -0.1819071952589948,
    -0.1821860859068651,
    -0.18231849212315632,
    -0.1823242188610244,
    -0.18232778855376458,
    -0.1823414279734353,
    -0.182981669840417,
    -0.1836121010677488,
    -0.18507077851171946,
    -0.18506509044202346,
    -0.18506412434687775,
    -0.185075151297879,
    -0.18506745509372877,
    -0.1856354165495409,
    -0.18425783615429422,
    -0.1829299052688484,
    -0.1825292314220611,
    -0.18252887876442728,
    -0.18224456052773833,
    -0.18171259062494008,
    -0.18186631959264415,
    -0.18200222221020865,
    -0.18169394528074242,
    -0.18194853358721796,
    -0.18248966881748965,
    -0.18064298916019067,
    -0.17987428825311386,
    -0.17979567367096563,
    -0.1801663864742297,
    -0.18015041323741451,
    -0.18027331918201892,
    -0.18032748603233537,
    -0.18046629871999065,
    -0.18066188460449684,
    -0.18077546779546866,
    -0.1813011646950154,
    -0.18122661226458656,
    -0.18113917716839842,
    -0.18123927454610406,
    -0.18111673225023253,
    -0.18125985757685242,
    -0.18135148224036282,
    -0.18125180858320455,
    -0.18132098504007899,
    -0.18128361401422421,
    -0.18108361224708092,
    -0.18063441104988262,
    -0.1811475538601168,
    -0.1811234359275665,
    -0.18113315406372074,
    -0.18115105285035848,
    -0.1810343271286979,
    -0.18115844656056285,
    -0.1809066784918677,
    -0.18113299416721737,
    -0.18021789689190285,
    -0.18066528643964666,
    -0.18056120444115625,
    -0.18053715766228212,
    -0.18062653221947392,
    -0.1808815583024251,
    -0.18116943363661075,
    -0.1814199232011651,
    -0.18169377951404972,
    -0.18157270527698932,
    -0.18213262307017242,
    -0.18127532066178412,
    -0.18103022229006174,
    -0.18034858132177592,
    -0.18045539555266504,
    -0.18001911184672095,
    -0.18023688906692567,
    -0.18087630156675447,
    -0.18076605658538916,
    -0.18127745354196997,
    -0.18146028182623555,
    -0.18134231315395438,
    -0.18122010987401463,
    -0.1817663674184319,
    -0.182547641106034,
    -0.18318912490504496,
    -0.18314908254206014,
    -0.18339165114842046,
    -0.18404933118273553,
    -0.18400259016932155,
    -0.18239419693059658,
    -0.18225477740312068,
    -0.18277845406261542,
    -0.18305911672009362,
    -0.1833079369646504,
    -0.18334664491739341,
    -0.18355254865824178,
    -0.18310378178610814,
    -0.18402131995468413,
    -0.18356306425785443,
    -0.1834816211279488,
    -0.18351379577771043,
    -0.18359086198639754,
    -0.18366599198632033,
    -0.1838523838857975,
    -0.18376746311901837,
    -0.18401453988966127,
    -0.18475075666691143,
    -0.18458284103708725,
    -0.18518732746303232,
    -0.1855905525731813,
    -0.18569101589819811,
    -0.18560713364162806,
    -0.18510743433842847,
    -0.18544506571474134,
    -0.18563152032254826,
    -0.185547488453682,
    -0.18550773747735222,
    -0.18541108322589173,
    -0.18570175485610502,
    -0.18502441103755426,
    -0.18627838621606246,
    -0.18747224820874503,
    -0.18733967559733727,
    -0.18741866255832584,
    -0.18746364210581734,
    -0.18769872841442012,
    -0.18863217725993248,
    -0.18690774753020722,
    -0.18728248197900207,
    -0.18850654138244632,
    -0.19007120715953624,
    -0.1901809145692643,
    -0.1893875136890912,
    -0.1891596803083007,
    -0.18938329622665792,
    -0.18993353681367103,
    -0.19010318374691004,
    -0.18980962527307377,
    -0.18996580897141738,
    -0.18985381918238692,
    -0.1892309768934635,
    -0.1895059247842441,
    -0.18833503377576738,
    -0.19043121508891314,
    -0.18591250293015607,
    -0.18559480306970974,
    -0.18872773348304364,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
]

PAM_EXPECTED_DISTANCES = [
    0.02313,
    0.02986,
    0.11968,
    0.12661,
    0.12679,
    0.12694,
    0.12743,
    0.13415,
    0.14256,
    0.29324,
    0.33685,
    0.34601,
    0.34643,
    0.34677,
    0.35077,
    0.35097,
    0.35112,
    0.35116,
    0.35141,
    0.36071,
    0.38831,
    0.40143,
    0.41098,
    0.42414,
    0.4627,
    0.46951,
    0.47431,
    0.47481,
    0.52863,
    0.53975,
    0.54158,
    0.56126,
    0.56677,
    0.80638,
    0.80791,
    0.81092,
    0.82953,
    0.83157,
    0.835,
    0.83984,
    0.84454,
    0.84501,
    0.84815,
    0.87542,
    0.88115,
    0.88122,
    0.89159,
    0.89234,
    0.91459,
    0.92201,
    0.93246,
    0.93313,
    0.93331,
    0.94031,
    0.9441,
    0.9442,
    0.9468,
    0.95163,
    0.9644,
    0.96519,
    0.96829,
    0.97543,
    0.97746,
    0.97881,
    0.98171,
    0.98213,
    0.98313,
    0.9908,
    0.99163,
    0.99878,
    1.00101,
    1.00392,
    1.00419,
    1.0098,
    1.01655,
    1.01783,
    1.02137,
    1.02173,
    1.02267,
    1.0235,
    1.02354,
    1.02395,
    1.03063,
    1.03199,
    1.03313,
    1.0332,
    1.04037,
    1.04223,
    1.04658,
    1.04759,
    1.04794,
    1.0503,
    1.05317,
    1.05393,
    1.05596,
    1.05603,
    1.0561,
    1.05793,
    1.06021,
    1.06078,
    1.06162,
    1.06241,
    1.06426,
    1.07221,
    1.07269,
    1.07482,
    1.08276,
    1.08421,
    1.0856,
    1.08691,
    1.0885,
    1.08945,
    1.09221,
    1.09233,
    1.09403,
    1.09788,
    1.10263,
    1.10407,
    1.10669,
    1.10976,
    1.11008,
    1.11071,
    1.11114,
    1.1186,
    1.11861,
    1.12023,
    1.12161,
    1.12242,
    1.12448,
    1.1245,
    1.13244,
    1.13527,
    1.14334,
    1.14733,
    1.14776,
    1.14805,
    1.15084,
    1.15249,
    1.15557,
    1.16235,
    1.16515,
    1.16525,
    1.16605,
    1.16731,
    1.16991,
    1.1703,
    1.17233,
    1.17255,
    1.17526,
    1.18319,
    1.18891,
    1.19525,
    1.20004,
    1.20038,
    1.20213,
    1.20267,
    1.20314,
    1.20526,
    1.2068,
    1.20739,
    1.21309,
    1.21563,
    1.21567,
    1.21753,
    1.21927,
    1.22026,
    1.22102,
    1.22103,
    1.22186,
    1.22477,
    1.22764,
    1.23767,
    1.2378,
    1.23786,
    1.23799,
    1.23821,
    1.23891,
    1.24064,
    1.2409,
    1.24265,
    1.24815,
    1.24824,
    1.24898,
    1.25012,
    1.25231,
    1.2552,
    1.25551,
    1.26021,
    1.26102,
    1.26119,
    1.26728,
    1.26818,
    1.27378,
    1.27386,
    1.27515,
    1.27745,
    1.28039,
    1.28342,
    1.28724,
    1.28928,
    1.29174,
    1.29628,
    1.29882,
    1.3078,
    1.31547,
    1.31739,
    1.3213,
    1.32264,
    1.3284,
    1.32901,
    1.34249,
    1.34664,
    1.34882,
    1.35137,
    1.35184,
    1.35393,
    1.35469,
    1.35798,
    1.36589,
    1.36858,
    1.3694,
    1.3695,
    1.37035,
    1.3707,
    1.37265,
    1.37327,
    1.37661,
    1.3824,
    1.38354,
    1.38994,
    1.39013,
    1.3912,
    1.39202,
    1.40245,
    1.40278,
    1.40427,
    1.40465,
    1.40472,
    1.4068,
    1.40814,
    1.42182,
    1.43825,
    1.44078,
    1.44111,
    1.44171,
    1.44221,
    1.44448,
    1.45444,
    1.46193,
    1.46377,
    1.48052,
    1.48169,
    1.48268,
    1.48998,
    1.49302,
    1.49656,
    1.50052,
    1.50344,
    1.50366,
    1.50558,
    1.50715,
    1.51227,
    1.51864,
    1.52501,
    1.57359,
    1.57728,
    1.61269,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
]


PAM_EXPECTED_SEQTODIS = [
    0,
    1,
    2,
    3,
    4,
    5,
    8,
    6,
    7,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    18,
    17,
    19,
    20,
    24,
    21,
    22,
    23,
    25,
    26,
    27,
    29,
    28,
    30,
    31,
    32,
    34,
    33,
    37,
    35,
    40,
    36,
    38,
    39,
    41,
    42,
    53,
    44,
    43,
    49,
    48,
    47,
    67,
    50,
    64,
    57,
    58,
    52,
    62,
    61,
    46,
    56,
    45,
    55,
    80,
    54,
    71,
    51,
    94,
    133,
    60,
    78,
    69,
    70,
    68,
    65,
    76,
    59,
    75,
    74,
    79,
    87,
    73,
    106,
    66,
    103,
    108,
    110,
    84,
    111,
    63,
    93,
    72,
    143,
    104,
    109,
    97,
    90,
    150,
    83,
    113,
    77,
    91,
    125,
    126,
    127,
    92,
    82,
    115,
    102,
    89,
    112,
    117,
    116,
    187,
    135,
    88,
    114,
    99,
    149,
    132,
    81,
    107,
    100,
    137,
    124,
    198,
    121,
    164,
    136,
    85,
    139,
    141,
    153,
    142,
    179,
    144,
    224,
    95,
    122,
    123,
    120,
    159,
    173,
    101,
    98,
    128,
    172,
    158,
    145,
    189,
    190,
    134,
    119,
    191,
    118,
    105,
    195,
    146,
    207,
    151,
    219,
    256,
    251,
    157,
    131,
    239,
    240,
    241,
    138,
    162,
    204,
    130,
    182,
    325,
    228,
    148,
    235,
    301,
    237,
    196,
    210,
    166,
    250,
    246,
    208,
    140,
    205,
    206,
    176,
    177,
    181,
    193,
    170,
    183,
    220,
    253,
    202,
    221,
    217,
    86,
    152,
    213,
    292,
    245,
    96,
    226,
    169,
    211,
    167,
    147,
    161,
    185,
    184,
    186,
    129,
    225,
    155,
    160,
    243,
    227,
    175,
    209,
    215,
    258,
    263,
    259,
    260,
    231,
    257,
    262,
    194,
    203,
    171,
    178,
    165,
    270,
    349,
    199,
    192,
    281,
    163,
    242,
    174,
    283,
    197,
    229,
    300,
    222,
    315,
    238,
    339,
    218,
    168,
    244,
    290,
    154,
    397,
    180,
    200,
    223,
    299,
    261,
    285,
    212,
    321,
    188,
    345,
    323,
    234,
    249,
    236,
    298,
    275,
    335,
    273,
    278,
    233,
    276,
    230,
    293,
    332,
    354,
    279,
    295,
    277,
    343,
    306,
    291,
    156,
    364,
    317,
    302,
    318,
    232,
    266,
    284,
    312,
    255,
    338,
    374,
    272,
    324,
    268,
    316,
    342,
    402,
    369,
    450,
    370,
    309,
    310,
    348,
    361,
    362,
    214,
    383,
    367,
    405,
    453,
    363,
    336,
    421,
    378,
    381,
    264,
    265,
    375,
    254,
    216,
    296,
    269,
    350,
    462,
    286,
    287,
    288,
    289,
    274,
    439,
    340,
    313,
    337,
    330,
    347,
    489,
    319,
    353,
    320,
    303,
    304,
    365,
    328,
    341,
    311,
    547,
    435,
    551,
    436,
    267,
    282,
    344,
    351,
    372,
    438,
    379,
    373,
    346,
    427,
    428,
    413,
    366,
    430,
    394,
    431,
    326,
    443,
    444,
    445,
    305,
    491,
    403,
    429,
    458,
    493,
    297,
    485,
    382,
    387,
    488,
    271,
    385,
    377,
    329,
    481,
    482,
    449,
    376,
    454,
    552,
    294,
    280,
    386,
    469,
    333,
    334,
    498,
    416,
    543,
    432,
    411,
    355,
    201,
    490,
    483,
    392,
    408,
    252,
    247,
    248,
    409,
    410,
    322,
    451,
    406,
    368,
    314,
    399,
    473,
    400,
    474,
    475,
    476,
    352,
    388,
    389,
    356,
    357,
    358,
    447,
    468,
    456,
    461,
    331,
    440,
    499,
    457,
    441,
    486,
    412,
    548,
    478,
    487,
    497,
    572,
    549,
    500,
    455,
    557,
    442,
    434,
    448,
    414,
    496,
    422,
    423,
    404,
    419,
    420,
    533,
    534,
    535,
    417,
    459,
    575,
    494,
    433,
    495,
    407,
    418,
    437,
    553,
    527,
    514,
    517,
    391,
    359,
    544,
    550,
    545,
    477,
    492,
    393,
    584,
    585,
    579,
    525,
    424,
    425,
    523,
    555,
    466,
    501,
    446,
    470,
    531,
    515,
    390,
    516,
    513,
    380,
    384,
    607,
    360,
    467,
    484,
    452,
    518,
    507,
    592,
    472,
    536,
    574,
    526,
    590,
    509,
    460,
    560,
    307,
    610,
    308,
    609,
    540,
    604,
    605,
    576,
    600,
    559,
    371,
    522,
    564,
    580,
    508,
    327,
    426,
    554,
    556,
    612,
    524,
    542,
    415,
    401,
    396,
    593,
    566,
    539,
    561,
    568,
    571,
    591,
    577,
    570,
    471,
    546,
    532,
    615,
    582,
    578,
    583,
    616,
    528,
    480,
    521,
    463,
    464,
    465,
    479,
    602,
    541,
    537,
    625,
    611,
    530,
    519,
    520,
    502,
    601,
    395,
    510,
    599,
    563,
    567,
    626,
    598,
    505,
    538,
    562,
    586,
    587,
    529,
    588,
    627,
    608,
    617,
    619,
    596,
    398,
    512,
    565,
    503,
    622,
    603,
    504,
    613,
    558,
    581,
    506,
    606,
    623,
    624,
    620,
    618,
    621,
    628,
    597,
    589,
    511,
    569,
    594,
    573,
    595,
    614,
]
