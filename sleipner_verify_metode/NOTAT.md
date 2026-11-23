# Kontroll av Sleipner-testen: er den rettferdig mot begge forklaringene?

Uavhengig kontroll, 19. november 2026 kveld. Alt ligger i `/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_verify_metode/`. Bare åpne data og eksperimentets egne filer er brukt. Eksperimentets filer er ikke endret.

## Dom

Regnestykkene i eksperimentet holder: jeg får de samme tallene fra stakkene, og min egen Zoeppritz-modell gir samme svar for det ene laget de testet. Men to av de fire konklusjonene er skrevet sterkere enn testen bærer:

- **Konklusjon 1 må være «kan gi», ikke «gir».** Underskuddet fra multiplene kommer av at multipler og primærer legger seg oppå hverandre i sonen, og det tåler bare noen få millisekunder forskyvning.
- **Konklusjon 2 gjelder bare ett homogent lag.** Lav Vs som følger de myke lagene i sonen, gir et underskudd av samme størrelse som det observerte.
- Konklusjon 3 holder. Konklusjon 4 holder, og er alvorligere enn rapporten sier.

Testen skiller altså ikke mellom (a) lav Vs og (b) multipler. Samlet dom: holder med rettelser.

## 1. Lav Vs: finnes det en geometri som gir underskuddet?

Egen modell: Zoeppritz i hver 1 m-grense (løst som 4×4-ligningssystem), egne stråler, gjennomsnitt over offsetene i hver stakk (mute 48°), 25 Hz Ricker. Sone/referanse for far/near, RMS-mål først og reflektormål etter skråstreken.

| Geometri | Vs −7 % | Vs −15 % | Vs −20 % |
|---|---:|---:|---:|
| Ingen endring (min modell) | 0,97 / 0,95 | | |
| Eksperimentets lag, 1358–1390 m | 0,96 / 0,90 | 0,99 / 0,90 | 1,03 / 0,91 |
| 30 geometrier med ett sammenhengende intervall (eksperimentets lag, 60, 90 og 150 m, hele sonen, gradvis topp eller bunn, over laget, under laget), RMS-mål | 0,91–1,14 | 0,87–1,34 | 0,86–1,45 |
| Alle bokser, topp 1200–1500 m, tykkelse 8–500 m (1064 stk.): laveste | 0,80 / 0,68 | | 0,65 / 0,56 |
| **Lav Vs i de myke lagene i sonen** (1353–1362, 1366–1373, 1375–1399 og 1425–1435 m) | **0,70 / 0,70** | **0,44 / 0,45** | **0,35 / 0,35** |
| Det samme i seks oppsett (25 Hz, 25 Hz med spøkelser og 35 Hz, hver med og uten NMO-strekk), RMS-mål | 0,62–0,75 | 0,36–0,61 | 0,26–0,58 |
| Ni måter å avgrense «de myke lagene» på (vindu 31, 61 eller 101 m, sonen eller ±30 m, gradert, lav Vp i stedet for lav impedans), RMS-mål | 0,70–0,81 | 0,44–0,70 | 0,35–0,68 |
| Det samme, mid/near, RMS-mål | 0,80–0,88 | 0,65–0,75 | 0,56–0,70 |

- Ett homogent lag gir ikke underskuddet. Det bekrefter eksperimentet for den geometrien de testet.
- «De myke lagene» er cellene der akustisk impedans ligger under glidende middel over 61 m. Det er leirpakken 1353–1399 m uten to tynne harde lag, pluss ett mykt lag 30 m dypere. Da faller far/near til 0,62–0,75 allerede ved 7 %, som er det de nærmeste skjærloggene tillater.
- Bare det å flytte bunnen av eksperimentets lag fra 1390 til 1396 m endrer svaret fra 0,99 til 0,87 (15 %, 25 Hz), og fra 0,89 til 0,74 med spøkelser og NMO.
- Denne modellen gir mindre underskudd på mid enn multippelmodellen gjør. Det ligger nærmere de ekte dataene.
- Forbehold: bare primærer, ikke FD. Med 18 Hz puls gir også eksperimentets lag 0,63–0,86.

## 2. Multiplene: hvor mye henger på antakelser?

Multiplene alene er FS − G i eksperimentets stakker. Jeg skalerer dem per stakk og måler med eksperimentets egne mål.

| Multipler | far/near sone/referanse | mid/near sone/referanse |
|---|---:|---:|
| Som i eksperimentet | 0,22–0,52 | 0,30–0,43 |
| Skalert til kopistyrken i de ekte dataene (0,73 / 0,62 / 0,55) | 0,57–0,62 | 0,45–0,52 |
| Skalert, og forskjøvet til forsinkelsene i de ekte dataene (−1 / −6 / −12 ms) | 0,69–0,80 | 0,64–0,75 |
| Som i eksperimentet, alle stakker forskjøvet −4 ms | 0,35–0,37 | 0,55–0,58 |
| Som i eksperimentet, alle stakker forskjøvet +4 ms | 0,65–1,04 | 0,37–0,52 |
| Som i eksperimentet, forskjøvet +10 ms | 3,6–4,8 | 2,4–2,9 |
| Observert | 0,48 | 0,94–1,00 eller 0,6–0,75 |

- Et underskudd står igjen etter skalering til observert styrke, men det er 0,6, ikke 0,48. Over de ni variantene med 25 Hz (mute, NMO, havbunn-R) blir det 0,27–0,63, og med 35 Hz 0,65–0,82.
- Multiplene virker mest ved å løfte near i sonen: near sone/referanse går fra 1,15 til 1,46, og far fra 1,19 til 0,91 (skalert). Lav Vs i de myke lagene lar near stå (1,14 til 1,08) og senker far (1,11 til 0,48 ved 15 %). Near og far hver for seg i de ekte dataene kan derfor skille de to forklaringene.
- Svaret avhenger av noen få millisekunder. Modellen binder til den ekte trasen med korrelasjon 0,45–0,59 og 1–6 ms forskyvning, så dette er innenfor usikkerheten.
- Glidende vindu på 84 ms gjennom 1100–1900 ms (66 vinduer, RMS-mål): uten havflate varierer vindu/referanse fra 0,51 til 1,56, og 9 % av vinduene ligger under 0,65. Med havflate 0,37–1,59, og 27 % under 0,65. Sonen er nummer 9 nedenfra, og det dypeste underskuddet ligger i vinduet 1360–1444 ms.
- Kopimålet (det filteret og sammenligningen med de ekte dataene bygger på) finner «første kopi» på −0,46 til −1,16 på near også i stakkene uten havflate, ved 86–142 ms. Styrken på kopiene er altså ikke et pålitelig mål. Forsinkelsen (110 ms og 220 ms) er det som peker på multipler i de ekte dataene.
- Havbunn-R 0,37–0,43 betyr lite (eksperimentets egne tall, bekreftet).

## 3. Mid

Svakere kopier på mid forklarer ikke avviket. Med skalerte multipler er mid/near 0,45–0,52. Av 2197 kombinasjoner av styrke per stakk (0 til 1,2) gir ingen både far/near 0,40–0,56 og mid/near over 0,90. Høyeste mid/near sammen med et slikt far-underskudd er 0,76–0,86, og da må multiplene være sterkest på far. Bare med multiplene forskjøvet −10 ms blir mid/near normal på RMS-målet (0,97–1,03; reflektormålet gir 0,73–0,79), og da er far/near 0,65–0,83.

## 4. Rapporten

Tallene i tabellene stemmer med `results.json`. Dette må rettes:

- `DX1_PLASSHOLDER` står igjen. 1 m-kjøringen finnes ikke i `results.json`.
- «0,14» i raden «reflektorer plukket på stakken med multipler» er medianen av to verdier. Den tredje reflektoren ga negativt forhold og ble forkastet.
- «Uten havflate er det aldri noe underskudd (0,94–1,26), uansett Vs»: med Vs-tilfellene er området 0,85–1,32.
- «og ikke nedover»: Vs −7 % gir 1,04 mot 1,07, altså litt ned.
- «alt fra 0,5 til 1,35» skal være 0,57.
- Courant-tallet i hovedkjøringene er 0,46, ikke 0,48.
- «Styrke 0,93 av Zoeppritz», «1,10 til 0,68» og «37 %» står ikke i `results.json` eller `tables.md`. Jeg får 0,995, 37 % og 1,13 til 0,70. 0,93 har jeg ikke gjenskapt.
- Kopitabellen sammenligner modellens kopier under én reflektor, målt på multiplene alene, med de ekte dataenes kopier målt på stakken selv. Målt på samme måte er modellen −0,94 / −0,52 / −0,75.
- «Morsøkten», «hovedmetoden», «den andre metoden», «FS» og «G» må forklares første gang.

## Slik er kontrollen sjekket

- Min Zoeppritz mot eksperimentets lukkede formel: største avvik 1,6 · 10⁻¹⁵ (5000 tilfeller). Energien er bevart til 1,2 · 10⁻¹⁵. Avvik fra Aki-Richards ved 2 % kontrast: 8 · 10⁻⁵.
- Min NMO-versjon med spøkelser mot eksperimentets Zoeppritz-stakker: korrelasjon 0,998 (near), 0,999 (mid), 0,982 (far).
- Første versjon hadde en feil: muten fjernet enkeltreflektorer nær 48° i NMO-varianten (far-korrelasjon 0,58). Den er rettet, og de gamle resultatene ligger i `v1_raw_angle_mute/`.
- Forskyvningene i punkt 2 er en grov etterligning: multippeltrasen flyttes i tid. Det er ikke en ny modellering.
- `p2`, `p3`, `p3b` og `p5` ble kjørt på den leide maskinen, resten på Linux-PC-en. `slp.py`, `model_1d.npz`, `stacks.npz`, `results.json` og `validation_images.json` er kopier av eksperimentets filer.

## Filer

| Fil | Innhold |
|---|---|
| `vz.py`, `selftest2.json` | Zoeppritz, stråler, stakker, mål, og selvtesten |
| `p1_lowvs.py`, `p1_results.json`, `p1.log` | Lav-Vs-geometriene, boksskann, fri profil |
| `p6_wavelet.py`, `p6_results.json` | Puls og NMO-strekk, sjekk mot eksperimentets stakker, andre Vs-lover |
| `p7_boundaries.py`, `p7_results.json` | Hvor topp og bunn av laget settes |
| `p2_scale.py`, `p2_results.json` | Skalerte multipler, glidende vindu |
| `p3_mid.py`, `p3b_grid_targets.py`, `p3_results.json`, `p3b_results.json`, `p3_grid3d.npy` | Mid, kopimålet uten multipler, 2197 styrkekombinasjoner |
| `p8_lags.py`, `p8_results.json` | Multipler med de ekte forsinkelsene, forskyvninger |
| `p4_audit.py`, `p4_audit.json`, `p5_checks.py`, `p5_checks.json` | Rapportens tall mot `results.json` |
