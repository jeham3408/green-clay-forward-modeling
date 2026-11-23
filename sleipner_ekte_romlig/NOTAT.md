# Sleipner: hvordan far/near-underskuddet varierer i rommet

Analyse «romlig» av de ekte stakkene (94p07 og 94p10), kjørt 19. november 2026. Alle tall står i `results.json` i samme mappe.

## Svaret først

**Den romlige oppførselen peker mot multipler (b), moderat sterkt. Den gir ingen støtte for leire med lav Vs (a), men utelukker ikke at 20–30 % av underskuddet er noe annet.**

- **Underskuddet følger kopien.** Blokk for blokk er far/near i sonen lavest der kopien på near er sterkest. Korrelasjonen er −0,48 i 2007-versjonen og −0,56 i 2010-versjonen, og usikkerheten når ikke null.
- **Det skjer fordi near går opp, ikke fordi far går ned.** Der kopien er sterk, er near i sonen 1,5 ganger sterkere enn der den er svak. Mid følger med (1,35 ganger). Far er nesten lik (0,92–0,93 ganger).
- **Der kopien er svakest, er underskuddet lite.** I tidelen av området med svakest kopi er far/near i sonen 0,76–0,82 av referansen når blokkene velges på uavhengige data, mot 0,54–0,57 i hele området. Kopien er ikke borte der (styrke 0,4), så resten kan ikke fordeles.
- **Test 1 (geometri) har ingen kraft.** Reflektorene i sonen ligger parallelt med de sterke reflektorene over (avstanden varierer med 2,5 ms standardavvik).
- **Test 3 (near mot far) er ikke avgjørende alene.** Mellom blokkene med dypest og grunnest underskudd er near 1,3–1,4 ganger så sterk og far 0,70 ganger så sterk, altså omtrent halvt om halvt.
- **Mid/near er lav i sonen** med en måling som ikke trenger innretting: 0,67–0,68 av referansen i begge prosesseringene. Det støtter den andre av de to tidligere målemetodene, ikke «normal mid/near». Far korrelerer nesten ikke med near i sonen (0,06–0,14 uten forskyvning, mot 0,71–0,73 ved reflektoren på 1316 ms), så en metode som bruker far som mellomledd, står svakt akkurat der.

## Hovedtall

Sone: 1400–1500 ms. Referanse: 1150–1360 ms og 1580–1880 ms. «Sone/referanse» er forholdet i sonen delt på forholdet i referansen. Tall i parentes er 95 % usikkerhet fra bootstrap over felt på 800 × 900 m.

| Måling | 94p07 | 94p10 |
|---|---|---|
| Blokker brukt (200 × 225 m), areal | 382 av 390, 17,2 km² | 381 av 390, 17,1 km² |
| far/near sone/referanse, hele området | 0,54 (0,52–0,58) | 0,54 (0,51–0,58) |
| Samme, bare referansen under sonen | 0,57 (0,54–0,61) | 0,59 (0,55–0,62) |
| mid/near sone/referanse | 0,67 (0,65–0,70) | 0,68 (0,66–0,71) |
| Kopi på near, mid, far (koeffisient ved forsinkelse) | −0,74 ved 110 ms, −0,38 ved 112, −0,26 ved 120 | −0,76 ved 110, −0,39 ved 112, −0,27 ved 120 |
| **Test 2:** korrelasjon underskudd mot kopistyrke | −0,48 (−0,63 til −0,40) | −0,56 (−0,64 til −0,45) |
| Samme, rangkorrelasjon | −0,60 (−0,66 til −0,51) | −0,60 (−0,67 til −0,50) |
| Samme, bare blokker minst 600 m fra hverandre (9 utvalg på ca. 42 blokker) | −0,55 (fra −0,72 til −0,43) | −0,58 (fra −0,68 til −0,45) |
| far/near sone/referanse: svak tredel / sterk tredel av kopien | 0,72 / 0,42 | 0,74 / 0,43 |
| Samme, valgt på uavhengig halvdel | 0,69 / 0,46 | 0,70 / 0,47 |
| Svakeste tidel, valgt på samme data | 0,96 (0,85–1,00) | 0,95 (0,78–1,08) |
| Svakeste tidel, valgt på uavhengig halvdel | 0,79 (0,73–0,87) | 0,76 (0,69–0,84) |
| Ekstrapolert til null kopi (samme data / uavhengig halvdel) | 0,93 (0,75–1,17) / 0,96 (0,88–1,05) | 0,85 (0,65–1,21) / 0,96 (0,87–1,06) |
| Nivå i sonen (mot referansen under), svak → sterk tredel: near | 0,83 → 1,21 | 0,82 → 1,21 |
| mid | 0,63 → 0,85 | 0,64 → 0,86 |
| far | 0,61 → 0,56 | 0,64 → 0,59 |
| Korrelasjon nivå mot kopistyrke: near / mid / far | +0,53 / +0,37 / −0,07 | +0,55 / +0,31 / −0,10 |
| far/mid sone/referanse (near er ikke med): alle / svakeste tidel på uavhengig halvdel | 0,81 / 0,96 (0,90–1,05) | 0,79 / 0,93 (0,86–1,01) |
| **Test 3:** dypeste mot grunneste fjerdedel, near / mid / far | 1,36 / 1,14 / 0,70 | 1,31 / 1,10 / 0,70 |
| Andel av variasjonen fra lav far / høy near | 0,51 / 0,49 | 0,54 / 0,46 |
| **Test 1:** avstand 1316 → 1386 ms, standardavvik og 5–95 % | 2,5 ms og 8,3 ms | 2,4 ms og 7,6 ms |
| Stigningstall (0 = følger de sterke, 1 = følger sonen) | 0,28 (−0,37 til 0,77) | −0,06 (−1,88 til 0,67) |
| Far forsinket etter near: ved 1277–1377 ms / 1407–1477 ms / 1697–1817 ms | 0,6–1,0 / 5,9–6,3 / 0,1–2,6 ms | 0,7–1,3 / 6,3–6,9 / 0,1–3,1 ms |

## Hva som ble gjort

1. **Blokkstakker.** Alle fire stakkene i begge prosesseringene ble lest med minnekart, åtte inlinjer om gangen, og midlet i celler på 8 × 9 traser (100 × 112,5 m). Fire celler utgjør en blokk på 200 × 225 m. Hele kuben er med: 3,1 × 5,85 km rundt 15/9-13 (249 × 468 traser med 12,5 m avstand). Seismikken ble aldri flyttet fra Linux-PC-en.
2. **Måling uten støyskjevhet.** Hver blokk deles i to halvdeler (cellene på diagonalen mot de to andre). Energien i et vindu er summen av produktet av de to halvdelene, så tilfeldig støy blåser den ikke opp. Amplitudeforholdet er roten av energiforholdet, og stakkene trenger ikke rettes inn mot hverandre.
3. **Reflektorer.** De sterke reflektorene ved 1276 og 1316 ms og sonens reflektorer ved 1386 og 1456 ms ble fulgt på fullstakken, blokk for blokk. Vinduene følger reflektoren ved 1316 ms (over) og 1386 ms (sonen og under).
4. **Kopien.** For hver blokk: minste kvadraters koeffisient mellom near i 1256–1356 ms og near 92–134 ms senere. For hele området gir det −0,74 ved 110 ms, som den tidligere analysen fant (−0,75).
5. **Testene.** Test 2 er korrelasjonen mellom underskuddet og kopistyrken, og nivået i grupper etter kopistyrke. Test 3 sammenligner hver stakk for seg mellom blokkene med dypest og grunnest underskudd, der rangering og måling er gjort på hver sin halvdel. Test 1 ser om profilen av far/near mot tid flytter seg med avstanden mellom reflektorene.

## Hva testene viser

**Test 1 har ingen kraft.** Avstanden mellom reflektoren ved 1316 ms og sonens reflektor ved 1386 ms varierer bare med 2,5 ms (standardavvik), mens strukturen selv har 30 ms relieff over området. Reflektoren ved 1386 ms ligger dessuten 110,9 ms (2007) og 111,2 ms (2010) under den ved 1276 ms, med motsatt fortegn. Det er nøyaktig kopiens forsinkelse, så det som ser ut som sonens egen reflektor på fullstakken, kan være kopien selv. Plukket ved 1456 ms hoppet mellom to hendelser (to topper rundt 20–25 ms fra hverandre i fordelingen) og er ikke brukt.

**Test 2 peker mot multipler.** Underskuddet er dypere der kopien er sterkere, i begge prosesseringene og også når bare blokker langt fra hverandre telles. Mekanismen er tydelig: near og mid i sonen stiger med kopistyrken, far gjør det ikke. Mid er egne data, så dette er ikke bare en regneeffekt på near.

**Test 3 er delt.** Halvparten av variasjonen i underskuddet mellom blokker kommer fra far. Men far varierer like mye mellom blokker i vinduer utenfor sonen (standardavvik av logaritmen 0,31–0,41 mot 0,32–0,35 i sonen), og den variasjonen henger ikke sammen med kopien. Den sier derfor lite om leire.

**Tillegg: far er forsinket i sonen.** Far ligger 6–7 ms etter near mellom 1407 og 1477 ms, mot rundt 1 ms ved de sterke reflektorene og 0–3 ms fra 1700 ms og ned. En multippel kommer senere på far enn på near (120 mot 110 ms). En vanlig refleksjon fra leire gjør ikke det. Dette er et indisium, ikke et bevis.

**Middelnivåene skiller ikke forklaringene.** Med modellens vinduer (1408–1492 ms mot 1100–1380 og 1560–1900 ms) måler vi far/near 0,54–0,55 og mid/near 0,66–0,67. Modellen (`prediksjoner.json`) gir 0,62 og 0,53 for skalerte multipler, og 0,56 og 0,70 for 15 % lavere Vs. Med offsetområdene som er tilpasset forsinkelsene, gir den 0,55 og 0,84 for multipler, og 0,57 og 0,97 for lav Vs. Ingen av dem treffer begge tallene, og ingen peker seg ut. Nivået for hver stakk alene er rundt halvparten av modellens (near 0,57–0,58 mot 1,1–1,5), så modellens «near alene» kan ikke brukes direkte. Det er den romlige oppførselen som skiller.

## Hva som er usikkert

- **Kopistyrken og far/near bruker de samme near-dataene.** Å velge blokker med svak kopi velger også blokker der near tilfeldigvis er svak, og det løfter far/near. Derfor faller svakeste tidel fra 0,95–0,96 til 0,76–0,79 når utvalget gjøres på en uavhengig halvdel. Det forsiktige tallet er 0,8.
- **Variasjonen i kopistyrke er trolig mest interferens** mellom kopien og lokale refleksjoner, ikke at multippelen selv er sterkere noen steder. Konklusjonen er den samme: amplituden på near i sonen styres av kopien.
- **2010-versjonen er ikke en uavhengig prøve.** Det er de samme feltdataene. Kartene over underskuddet korrelerer 0,91 mellom prosesseringene.
- **Et restunderskudd på 20–30 % kan ikke utelukkes.** Ekstrapolert til null kopi er far/near 0,85–0,96 av referansen, med usikkerhet fra 0,65 til 1,21. Et underskudd på 40–45 % fra leire alene passer dårlig med svakeste tidel (0,76–0,82 selv med kopi til stede).
- **Stakkede data kan ikke skille multipler og primærrefleksjoner direkte.** Forsinkelsen på far er det nærmeste vi kommer.

## Hva som ville endret konklusjonen

- Gathers som viser at hendelsene i sonen flater ut som primærrefleksjoner.
- At underskuddet på rundt 20 % i blokkene med svakest kopi står igjen etter fjerning av multipler på gathers.
- En skjærlogg gjennom laget som viser Vs 10–15 % under trenden.

## Filer

- `s1_blockstack.py`: blokkstakker (`blockstack.npz`)
- `s2_analyse.py`: alle målinger og tester (`results.json`, `maps.npz`, `s2.log`)
- `s3_figs.py`: figurene
- `s5_derived.py`: legger de få avledede tallene i notatet inn i `results.json` (kjøres etter `s2_analyse.py`)
- `fig_romlig_kart.png`: kart for 2007-versjonen (tid, avstand mellom reflektorer, kopistyrke, far/near, near og far hver for seg)
- `fig_romlig_tester.png`: underskudd mot kopistyrke, forhold mot tid, og nivå per stakk etter kopistyrke (gruppene i søylediagrammet er valgt på samme data; spredningsplottet bruker referansen under sonen, der korrelasjonen er −0,42 og −0,44)
