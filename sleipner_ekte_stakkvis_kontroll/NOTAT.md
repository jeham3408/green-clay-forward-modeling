# Kontroll av «sleipner_ekte_stakkvis» (23.11, ca. 23.25)

## Svar først
**Dommen holder med rettelser.** Jeg regnet om de viktigste tallene med egen kode (egne 12 ms-binner, ikke deres funksjoner) og får det samme: near i sonen 1,40 av referansen, mid 1,07–1,09, far 0,94–0,95, mid/near 0,76–0,78, far/near 0,67, i begge prosesseringene. Avvik fra deres tall er høyst 0,01. Jeg fant ingen akse-, tidsnull- eller gainfeil som kan lage mønsteret.

Det som holder: **near er for sterk i sonen, ikke far for svak.** Near-forsterkningen sitter i ett vindu (1420–1456 ms, near 1,62 mot mid 1,11 og far 0,98), omtrent 108 ms under den sterke pakken ved 1312–1348 ms. Det er der en vannlagsmultippel skal ligge, og lav Vs kan ikke gjøre near sterkere. Når filteret E fjerner kopiene, faller near fra 1,40 til 1,12–1,13 mens far bare går fra 0,94 til 0,89.

Det som er for sterkt formulert hos dem:
1. «Mid/near er avgjort lav» er riktig for de stabile målene (RMS 0,76–0,78 og minste kvadrater 0,76–0,78 ufiltrert), men mid/near skiller ikke mellom forklaringene. Modellen gir 0,53–0,84 for multipler og 0,70–0,97 for 15 % lavere Vs, og målt 0,76 ligger i begge områdene. Det lar seg ikke bruke som bevis for (b).
2. Etter filteret er mid/near ikke avgjort: RMS 0,88–0,89, minste kvadrater 0,67–0,69, og instrumentmålet 0,94–1,00. Så «0,86–0,89 etter filter» er ett av tre sprikende mål.
3. Far/near er et dårlig mål her, ikke bare «0,47–0,72». Korrelasjonen near–far i sonen er 0,06–0,09 (0,29 i referansen), så alle mål som deler på near·far eller regresserer far på near blir meningsløse (minste kvadrater 0,15–0,19). Bare RMS-forholdet (0,67) er stabilt.
4. Selv uten Vs-effekt forklarer multipler ikke hele mønsteret like pent: modellen gir med skalerte multipler near 1,46, mid 0,77 eller 1,22 (avhengig av offsetvalg), far 0,81–0,91. Real far 0,94 er litt høyere enn begge. Det er bra for (b), men kombinasjonen 7 % lavere Vs pluss multipler gir far 0,61–0,66 i modellen, altså mye lavere enn målt 0,94. Mønsteret **trenger altså ikke** en Vs-effekt, og en stor Vs-effekt (15 %) på toppen av multipler er utelukket om modellens nivå er riktig.
5. Resten etter filter (far/near 0,79, mid/near 0,87–0,89) er 2,5–3 standardavvik fra andre vinduer (bare 38 overlappende vinduer, ca. 11 uavhengige), og sonen er plukket fordi den er lav. Den resten bør kalles «uforklart», ikke «ligner 7 % lavere Vs».

**Hvilken forklaring dataene peker mot:** multipler (b) som hovedårsak til near-overskuddet og dermed det meste av far/near-underskuddet. Styrke: moderat til sterk for at multipler finnes og løfter near; svak for at det ikke i tillegg ligger noe i resten. Lav Vs er ikke vist og ikke utelukket for en liten rest.

## Hovedtall (mine, begge prosesseringer; hele kuben 249×468 traser; delområder = 18 stk. à ca. 1 km²)
| mål (sone/referanse) | 94p07 | 94p10 | deres 94p07/94p10 |
|---|---|---|---|
| near RMS | 1,404 | 1,406 | 1,40 / 1,40 |
| near middel |amplitude| | 1,423 | 1,424 | (1,40–1,55 over mål) |
| near median over traser | 1,35 | 1,35 | |
| near delområder (median, SE, min–maks) | 1,41, 0,011, 1,29–1,48 | 1,41, 0,011, 1,31–1,48 | 1,41, 0,011 |
| mid RMS | 1,07 | 1,09 | 1,07 / 1,09 |
| far RMS | 0,938 | 0,946 | 0,94 / 0,95 |
| mid/near RMS | 0,762 | 0,776 | 0,77 / 0,78 |
| mid/near delområder | 0,76 ± 0,009 (0,71–0,85) | 0,77 ± 0,009 (0,72–0,86) | 0,76 ± 0,009 |
| mid/near minste kvadrater (Σnm/Σnn) | 0,756 | 0,782 | |
| mid/near instrumentmål (far som hjelp) | 1,65 | 2,19 | 1,68 / 2,50 |
| far/near RMS | 0,668 | 0,673 | 0,67 / 0,68 |
| far/near delområder | 0,66 ± 0,011 | 0,68 ± 0,010 | 0,67 ± 0,011 |
| far/near minste kvadrater (Σnf/Σnn) | 0,19 | 0,15 | 0,47–0,49 (deres LS per reflektor) |
| korrelasjon near–far, sone / referanse | 0,085 / 0,294 | 0,062 / 0,285 | 0,11 / 0,41 (blokktraser) |
| korrelasjon near–mid, sone / referanse | 0,46 / 0,46 | 0,45 / 0,45 | 0,63 / 0,61 |
| near RMS etter filter E | 1,12 | 1,13 | 1,13 |
| mid RMS etter filter E | 0,98 | 1,01 | 0,99 |
| far RMS etter filter E | 0,89 | 0,89 | 0,89 |
| mid/near etter E (RMS / LS / instrument) | 0,875 / 0,67 / 0,94 | 0,889 / 0,69 / 1,00 | 0,86–0,89 |
| far/near etter E (RMS / LS) | 0,79 / 0,49 | 0,79 / 0,47 | 0,77–0,79 |

(Mine er ratio av pooled RMS i 7 binner à 12 ms: sone 1408–1492 ms; referanse 1108–1372 og 1564–1900 ms. Forskjellen til deres vinduer er under 0,01.)

## Tidsprofil (36 ms-vinduer, RMS relativt til referansen, 94p07 ufiltrert)
Near / mid / far: 1384: 1,13/0,96/0,87; **1420: 1,62/1,11/0,98**; 1456: 1,21/1,06/0,91; 1492: 1,11/1,08/0,85; 1528: 1,14/1,06/0,86; 1564: 1,01/1,00/0,97. Den sterke pakken 1312–1348: 1,64/1,70/1,72 (alle tre like, så ingen gainforskjell der). Far/near går ned til 0,61 ved 1420 og er 0,75–0,76 ved 1456–1528. Etter filter E: far/near 0,77–0,84 i 1384–1528, mid/near 0,84–0,98. Underskuddet ligger altså **ikke** bare på leirlaget (1412–1446 ms) men bredt 1384–1530 ms (150 ms), og delvis i en annen vannlagskopi (+220 ms ≈ 1536 ms).

## Uavhengig kontroll av kopien
Regresjon av vindu 1366–1466 ms på samme stakk 110 ms tidligere (kun «vinduet over», ikke den spesifikke pakken, derfor lavere enn deres 49–52 %): near −0,54 ved 112 ms (forklarer 29 % av energien), mid −0,29 ved 114 ms (15 %), far −0,21 ved 120 ms (10 %). Samme rekkefølge near > mid > far og voksende forsinkelse med offset (moveout) som hos dem. Men en slik negativ kobling finnes også i andre vinduer (1500–1600: near −0,36, 20 %; 1600–1700: −0,31, 16 %; 1200–1300: −0,26, 6 %), så kopien er ikke spesifikk for sonen. Den er bare mest effektiv der den sterke pakken ligger rett over.

## Gain
RMS per 96 ms (94p07, near/mid/far): 1000 ms: 0,18/0,15/0,10; 1192: 0,21/0,19/0,13; 1288: 0,28/0,25/0,18; 1384: 0,28/0,19/0,12; 1480: 0,23/0,19/0,11; 1672: 0,18/0,15/0,11. Forholdene mid/near og far/near er jevne (0,85–0,88 og 0,60–0,65) utenfor 1384–1530 ms og faller bare der. Et glatt gain kan ikke lage en 150 ms bred dupp med rett kant. Gain forklarer ikke mønsteret. Tidsinnretting: beste forsinkelse near–mid 2 ms (1 sample) i sone og referanse; near–far 2 ms i referansen men 6 ms i sonen (der korrelasjonen er 0,06–0,09, så toppen er ikke meningsfull); mid–far 2–4 ms. Under én sample for stabile par.

## Forsøk på å knekke konklusjonen
- **Akse/tidsnull/brønn:** resultatene er for hele kuben og 18 delområder (alle har mid/near 0,71–0,86 og far/near 0,59–0,74); brønnposisjonen betyr dermed ingenting for hovedtallene. Ikke kontrollert: selve brønnkoblingen til kubens akser.
- **Utvalgseffekt:** sonen ble valgt fordi far/near var lav. Mot 38 overlappende sammenligningsvinduer i referansen ligger far/near under alle (placebo-snitt 1,00, sd 0,055–0,07), men siden vinduet er plukket, teller ikke dette som uavhengig. Mid/near (0,76) var *ikke* valgt (første analyse fant 0,94–1,00), og er også lavere enn alle placebo-vinduer (snitt 0,99, sd 0,023–0,028). Det er uavhengig bekreftelse på at noe skjer, men ikke på hva.
- **Absolutt near-nivå:** near alene varierer mye mellom vinduer i referansen (0,79–1,55, sd 0,21). 1,40 er bare 2 sd over snittet, og et enkelt vindu i referansen når 1,55. Near-nivået i seg selv er derfor *ikke* et bevis. Det som veier tungt er (i) at økningen sitter nøyaktig 108–112 ms under en sterk pakke og bare på near, (ii) at filtrene tar den bort, (iii) at mid og far ikke følger. Dette hviler ikke på modellens nullnivå, men på mønsteret.
- **RMS trekkes mot 1 av støy.** Delvis nevnt hos dem. Enig.
- **Referansevalg:** referansen inneholder den sterke pakken 1276–1348 (1,2–1,7) og stille sene vinduer (0,8–0,9). Det flytter alle tre stakkers absolutte nivåer likt og endrer ikke forholdene mellom stakkene.

## Hva som er usikkert
- Modellens nullnivå binder bare svakt til ekte data (korrelasjon 0,45–0,59), så sammenligning av 1,40 mot 1,15 (modell uten forklaring) er rettledende.
- Filtre på stakker rydder ikke pålitelig opp (vist av den syntetiske testen). Tallene etter filter er derfor nedre grenser for hva multipler gjør.
- 94p07 og 94p10 er ikke uavhengige målinger av samme ting (samme feltdata, toppkorrelasjon 0,87–0,94). At de gir samme tall er ingen replikasjon.
- Hvis jeg skulle skifte mening mot lav Vs: near i sonen *uten* kopier måtte vise seg tydelig lavere enn 1,12 samtidig som far/near forble rundt 0,8, eller et filter som er validert på syntetiske data måtte la far/near stå uendret mens near gikk ned. Mot multipler: kopistyrken i sonen under reflektorer som *ikke* har en sterk pakke 110 ms over, burde være like stor som under pakken (den er det ikke: kopien er dårligst i 1200–1300).

## Filer
Scripts: `kontroll.py` (binner), `ana.py`, `ana2.py`, `fig.py`. Resultater: `results.json`, `results_ekstra.json`. Figur: `fig_kontroll_profil.png`. Mappen inneholdt fra før filer jeg ikke har laget (`acc_*.npz`, `post*.py`, `post*.log`, `check0.*`, `fig_kontroll_kopi_og_sted.png`, `fig_kontroll_tid.png`); en tidligere `results.json` kan være overskrevet av min.
