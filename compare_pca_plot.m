
if ~exist('pca_single1')
    pca_single1=[0.00307543179951608180999755859375000000000000000000000,0.00263050035573542118072509765625000000000000000000000,0.00645024096593260765075683593750000000000000000000000,0.05646473541855812072753906250000000000000000000000000,0.02166110277175903320312500000000000000000000000000000,0.00773493433371186256408691406250000000000000000000000,0.00289984815753996372222900390625000000000000000000000,0.00360553315840661525726318359375000000000000000000000,0.00302362535148859024047851562500000000000000000000000,0.00724868010729551315307617187500000000000000000000000,0.00936621706932783126831054687500000000000000000000000,0.00704087642952799797058105468750000000000000000000000,0.00419323705136775970458984375000000000000000000000000,0.00603672908619046211242675781250000000000000000000000,0.00203026295639574527740478515625000000000000000000000,0.00334263034164905548095703125000000000000000000000000,0.00401877192780375480651855468750000000000000000000000,0.00481040123850107192993164062500000000000000000000000,0.01629068888723850250244140625000000000000000000000000,0.00496354745700955390930175781250000000000000000000000,0.00745899928733706474304199218750000000000000000000000,0.00266993977129459381103515625000000000000000000000000,0.00468170735985040664672851562500000000000000000000000,0.00456433650106191635131835937500000000000000000000000,0.01038312539458274841308593750000000000000000000000000,0.00356816477142274379730224609375000000000000000000000,0.00303993653506040573120117187500000000000000000000000,0.00287852878682315349578857421875000000000000000000000,0.00400106329470872879028320312500000000000000000000000,0.00274969055317342281341552734375000000000000000000000,0.00317430309951305389404296875000000000000000000000000,0.00629630777984857559204101562500000000000000000000000,0.00874250568449497222900390625000000000000000000000000,0.10941506922245025634765625000000000000000000000000000,0.00704983202740550041198730468750000000000000000000000,0.00627906760200858116149902343750000000000000000000000,0.00609097862616181373596191406250000000000000000000000,0.00433360645547509193420410156250000000000000000000000,0.00480462703853845596313476562500000000000000000000000];
    pca_single2=[-0.00000635570222584647126495838165283203125000000000000,0.00002838833825080655515193939208984375000000000000000,-0.00004350569361122325062751770019531250000000000000000,0.31456625461578369140625000000000000000000000000000000,-0.00061957107391208410263061523437500000000000000000000,-0.00010526354162720963358879089355468750000000000000000,0.00006960846803849563002586364746093750000000000000000,-0.00002191482781199738383293151855468750000000000000000,0.00001458738915971480309963226318359375000000000000000,0.00007043837103992700576782226562500000000000000000000,0.00023906696878839284181594848632812500000000000000000,0.00013605452841147780418395996093750000000000000000000,-0.00007072330481605604290962219238281250000000000000000,-0.00006837375258328393101692199707031250000000000000000,0.00000694067739459569565951824188232421875000000000000,-0.00006892944657010957598686218261718750000000000000000,-0.00006048723662388511002063751220703125000000000000000,-0.00014080690743867307901382446289062500000000000000000,-0.00018039403948932886123657226562500000000000000000000,-0.00010482132347533479332923889160156250000000000000000,0.00045473370118997991085052490234375000000000000000000,-0.00000765615823183907195925712585449218750000000000000,0.00003005178405146580189466476440429687500000000000000,-0.00006924525951035320758819580078125000000000000000000,-0.00031998855411075055599212646484375000000000000000000,0.00005944378426647745072841644287109375000000000000000,-0.00002885763751692138612270355224609375000000000000000,-0.00000863039986143121495842933654785156250000000000000,0.00003540999750839546322822570800781250000000000000000,-0.00001792297007341403514146804809570312500000000000000,0.00008225391502492129802703857421875000000000000000000,0.00004450967389857396483421325683593750000000000000000,-0.00006195052264956757426261901855468750000000000000000,-0.65451043844223022460937500000000000000000000000000000,-0.00006894439866300672292709350585937500000000000000000,-0.00002711018169065937399864196777343750000000000000000,0.00007996537169674411416053771972656250000000000000000,0.00001907948535517789423465728759765625000000000000000,-0.00013401798787526786327362060546875000000000000000000];
    pca_double1=[0.00187432133874235423535992151045093123684637248516083,0.00184287925483541249575825116124860869604162871837616,0.00203495207394839474485115360380405036266893148422241,0.02280231707120123837984060344297176925465464591979980,0.01541403346346682454170995413278433261439204216003418,0.00337060095995851292591249936947406240506097674369812,0.00197151718168999019742404499311305698938667774200439,0.00221118790024814803815700336997451813658699393272400,0.00207502517609811512958040147225347027415409684181213,0.00423718752556855818985903638917989155743271112442017,0.00533679995407691787345783041018876247107982635498047,0.00556854574132827854354710339634948468301445245742798,0.00174672442280091955515985979729975952068343758583069,0.00164916084178119812415908018721211192314513027667999,0.00125623794613190426817939115977651454159058630466461,0.00202021887682253369067364445754719781689345836639404,0.00212896979886713138760412356020879087736830115318298,0.00218932969656861512938728075994276878191158175468445,0.00715175497171101348659050245260004885494709014892578,0.00326507231025674037744321864806806843262165784835815,0.00319026173995537584357484739427945896750316023826599,0.00174095618195268392325303352663468103855848312377930,0.00170030419354487470975922924765200150432065129280090,0.00172593599949879953660281106664342587464489042758942,0.00515533566642469194551967603956654784269630908966064,0.00195045257744091288243382109612866770476102828979492,0.00192690710283112626513557508189933287212625145912170,0.00175495065625294838045078460453396473894827067852020,0.00195190288969025781476429592942167801083996891975403,0.00194235112464426857778465596027217543451115489006042,0.00178538830280340556251605921289637990412302315235138,0.00292202099327308400730784931909056467702612280845642,0.00353098968134116274544820335279382561566308140754700,0.05876867851592301739138690663821762427687644958496094,0.00395582855521717496188882634555739059578627347946167,0.00217495888834310854842435389855381799861788749694824,0.00333686489460178732982531535355974483536556363105774,0.00320143412752413444519183904901638015871867537498474,0.00317165206952416513577763801379205688135698437690735];
    pca_double2=[0.00000000001809300535527640459325986216974592934655930,0.00000000001747745697713128295006916550316148446439501,0.00000000002131936468428566421521476972283199672469900,0.00000000276906903024653399948201594330584457415156407,0.00000000123225285650262710104172436643369531461544852,0.00000000005838654572242841846665813543322954941305891,0.00000000002019595127917901512522172582114895166005564,0.00000000002547104231300585548210640276799344864158692,0.00000000002235813660273207855269971444151953440804270,0.00000000009367476666650535029475800372186951314734671,0.00000000014814100438732856547029561441119419745726660,0.00000000016189296292305608005040422664493537582841576,0.00000000001585550418625716641568430190844743313607901,0.00000000001412433657666547688501046296330826098665290,0.00000000000822213597622815281996548350794214621468925,0.00000000001916320533792738179035012136298683597945836,0.00000000002353355717324051668753430193274520576193765,0.00000000002498020937647652656713484676515967929794071,0.00000000026017010217000351559207558752994261103141937,0.00000000005628745743174293914630904610769447388857190,0.00000000005346397649362491917191233600346216611048922,0.00000000001558518667325031171629602546390916265398041,0.00000000001488288345507439769326278383543380335468642,0.00000000001533507241712195465874559955245378609772766,0.00000000013630026339585068171728234523291424518109949,0.00000000001952134621024381040487506561465364798499400,0.00000000001913699005453399388554333894531278049966350,0.00000000001589981819193208396868107567047095820406133,0.00000000001970523045523151054308490360249131803649369,0.00000000001952479684134558322502665666008044064241855,0.00000000001648135314575471042548206564981990102439213,0.00000000004407804318212711425166941858680826147459664,0.00000000006521711613951682959337304058437595130626185,0.00000002634172075823844234566997984502839669218587915,0.00000000008154137183204583698816696435114972918434262,0.00000000002415296743883607719119533538599018480594327,0.00000000005748700787732190144375947115803723082988475,0.00000000005269786379533034325218397821830496226269691,0.00000000005219566752459378693162063101147025382525735];
    mir_single=[2.05665211303340811355155892670154571533203125000000000,1.88125867688745529449079185724258422851562500000000000,2.73149622880600873031653463840484619140625000000000000,5.37732579232232410504366271197795867919921875000000000,7.74423655625497531218570657074451446533203125000000000,4.87177971975603441023849882185459136962890625000000000,4.19170452366563495161244645714759826660156250000000000,7.51558917058895303853205405175685882568359375000000000,1.22355544356480550050036981701850891113281250000000000,6.83910172157294482531142421066761016845703125000000000,6.20395828133939630788518115878105163574218750000000000,6.29960450753611667096265591681003570556640625000000000,3.28632970981209382443921640515327453613281250000000000,1.30357765298987260393914766609668731689453125000000000,1.21858705638362607714952901005744934082031250000000000,5.87647949837537453277036547660827636718750000000000000,3.65883807168825114786159247159957885742187500000000000,5.84823830306578429372166283428668975830078125000000000,9.06534658901460943525307811796665191650390625000000000,6.86065845386502815017593093216419219970703125000000000,6.47633780818756576991290785372257232666015625000000000,5.31109395815474272239953279495239257812500000000000000,2.71452863135226607482763938605785369873046875000000000,2.67810850342721096239984035491943359375000000000000000,4.83958845266215575975365936756134033203125000000000000,2.95762691151310264103813096880912780761718750000000000,2.25376786286722108343383297324180603027343750000000000,4.23244296722862145543331280350685119628906250000000000,5.31516490314413658779812976717948913574218750000000000,3.49584113630140791428857482969760894775390625000000000,6.27903674558729107957333326339721679687500000000000000,7.72222779404927450741524808108806610107421875000000000,12.61814863717944490417721681296825408935546875000000000,-13.47604823506856064341263845562934875488281250000000000,8.18175106809735552815254777669906616210937500000000000,1.45301950292207493475871160626411437988281250000000000,4.74674154540514336986234411597251892089843750000000000,7.90136033454058406277908943593502044677734375000000000,5.26413736391668862779624760150909423828125000000000000];
    mir_double=[2.05894692450459615429281257092952728271484375000000000,1.89120967456943844808847643435001373291015625000000000,2.72907822166808955444139428436756134033203125000000000,5.37247249080877509186393581330776214599609375000000000,7.74264502346761673834407702088356018066406250000000000,5.07899032493941149368765763938426971435546875000000000,4.19235829722276776010403409600257873535156250000000000,7.51147016580171111854724586009979248046875000000000000,1.22275450642570149284438230097293853759765625000000000,6.85283961617903969454346224665641784667968750000000000,6.22907681494757525797467678785324096679687500000000000,6.34443022725361061020521447062492370605468750000000000,3.48441372845070418406976386904716491699218750000000000,1.30614039168870021967450156807899475097656250000000000,1.21442257497045602576690725982189178466796875000000000,5.88414139584898521206923760473728179931640625000000000,3.64656350009039442738867364823818206787109375000000000,5.81388303442020060174399986863136291503906250000000000,9.06515463666977439061156474053859710693359375000000000,6.85608132705169737164396792650222778320312500000000000,6.56778997529045227565802633762359619140625000000000000,5.31186640436413881616317667067050933837890625000000000,2.71710666498529462842270731925964355468750000000000000,2.67587337355558929630205966532230377197265625000000000,4.82358450270015737260109744966030120849609375000000000,2.95287322680917441175552085041999816894531250000000000,2.25289461767545162729220464825630187988281250000000000,4.23388129171286209384561516344547271728515625000000000,5.30133567471739297616295516490936279296875000000000000,3.50210962135696490804548375308513641357421875000000000,6.26825717019858075218508020043373107910156250000000000,7.73340773416424553943215869367122650146484375000000000,12.61767176235878196166595444083213806152343750000000000,-13.31845596991922775487182661890983581542968750000000000,8.22697071631523613177705556154251098632812500000000000,1.45278401336531715060118585824966430664062500000000000,4.76688578565097031969344243407249450683593750000000000,8.06413763849667475369642488658428192138671875000000000,5.29039615769187321347999386489391326904296875000000000];
    mir_double2=[38.40615942969313323374080937355756759643554687500000000,40.38588175508974131844297517091035842895507812500000000,43.50138152194136864636675454676151275634765625000000000,49.34211331497627384123916272073984146118164062500000000,41.18803300079821383405942469835281372070312500000000000,24.44763079282367357336624991148710250854492187500000000,37.95082456105247104005684377625584602355957031250000000,39.00619873099279288908292073756456375122070312500000000,40.22151934609554757571459049358963966369628906250000000,54.63779501963047380286297993734478950500488281250000000,70.46930232544798400340368971228599548339843750000000000,72.97113477192746699984127189964056015014648437500000000,30.73357279726663193741842405870556831359863281250000000,32.71769379397355237415467854589223861694335937500000000,23.14536837194680174434324726462364196777343750000000000,44.12582222392654784925980493426322937011718750000000000,50.50766729129650656204830738715827465057373046875000000,50.72494977431631468789419159293174743652343750000000000,55.06728871390716051337221870198845863342285156250000000,41.22090731466917645775538403540849685668945312500000000,43.10390271129578820819006068632006645202636718750000000,55.83858591698236750744399614632129669189453125000000000,55.70212326590208107290891348384320735931396484375000000,56.70234311843567809319210937246680259704589843750000000,78.56709537083993666328751714900135993957519531250000000,38.43289025898043576034979196265339851379394531250000000,38.87103894058559205859637586399912834167480468750000000,31.38612864970468763203825801610946655273437500000000000,36.05001067856309759918076451867818832397460937500000000,42.08574299906081250810530036687850952148437500000000000,18.16409761840662895338027738034725189208984375000000000,40.85707861119510653225006535649299621582031250000000000,32.54690576742706298318807967007160186767578125000000000,227.96177593875526667943631764501333236694335937500000000,58.24345740513842883956385776400566101074218750000000000,51.10069365798104001896717818453907966613769531250000000,67.79438251959871308827132452279329299926757812500000000,66.40234351801504431023204233497381210327148437500000000,72.43450857577931856212671846151351928710937500000000000];
end

pca_single1 = abs(pca_single1);
pca_single2 = abs(pca_single2);
pca_double1 = abs(pca_double1);
pca_double2 = abs(pca_double2);

fprintf('\nLast eigenv using PCA\n')
fprintf('Mean and std single: %1.4f (%1.4f)\n', mean([pca_single1]), std([pca_single1]));
fprintf('Mean and std double: %1.4f (%1.4f)\n', mean([pca_double1]), std([pca_double1]));
fprintf('Sign test: %g\n', signtest([pca_single1]-[pca_double1]));

fprintf('\nLast eigenv using PCA of cov matrix\n')
fprintf('Mean and std single: %1.10f (%1.10f)\n', mean([pca_single2]), std([pca_single2]));
fprintf('Mean and std double: %1.10f (%1.10f)\n', mean([pca_double2]), std([pca_double2]));
fprintf('Sign test: %g\n', signtest([pca_single2]-[pca_double2]));

fprintf('\nMIR\n');
fprintf('Mean and std MIR single : %1.4f (%1.4f)\n', mean([mir_single]), std([mir_single]));
fprintf('Mean and std MIR double : %1.4f (%1.4f)\n', mean([mir_double]), std([mir_double]));
fprintf('Mean and std MIR double2: %1.4f (%1.4f)\n', mean([mir_double2]), std([mir_double2]));
fprintf('Mean and std MIR double3: %1.4f (%1.4f)\n', mean([mir_double3]), std([mir_double3]));

fprintf('\nPMI\n');
fprintf('Mean and std PMI single : %1.8f (%1.4f)\n', mean([pmi_single]), std([pmi_single]));
fprintf('Mean and std PMI double : %1.8f (%1.4f)\n', mean([pmi_double]), std([pmi_double]));
fprintf('Mean and std PMI double2: %1.8f (%1.4f)\n', mean([pmi_double2]), std([pmi_double2]));
fprintf('Mean and std PMI double3: %1.8f (%1.4f)\n', mean([pmi_double3]), std([pmi_double3]));


return



figure; 
hist2([pca_single1], [pca_double1]);
legend({'single', 'double'})

figure; 
plot([pca_single], [pca_double], 'b.');
[ypred, alpha, rsq, slope, intercept] = fastregress([pca_single], [pca_double], 1, 0);
axis equal
legend({'single', 'double'})