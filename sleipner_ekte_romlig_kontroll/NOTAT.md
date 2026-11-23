# Kontroll av "sleipner_ekte_romlig" (romlig analyse av far/near-underskuddet)

Alle tall under kommer fra egen kode (k1_cells.py, k2_analyse.py, k3_extra.py) og ligger i results.json i denne mappa. Koden bruker ikke funksjoner fra den andre analysen.

## Svar først

**Tallene holder, men konklusjonen er ordlagt sterkere enn dataene tillater.** Jeg får de samme hovedtallene med en annen metode, men tre av de fire testene skiller ikke forklaringene så godt som notatet gir inntrykk av.

- Det som holder: kopien 110 ms under de sterke reflektorene finnes overalt (koeffisient -0,72 til -0,76), far/near i sonen er 0,55 av referansen (samme som deres 0,54), mid/near 0,69 (deres 0,67). Near i sonen er omtrent 1,45 ganger referansen mens far er 0,8-0,9. Det ligner mye på modellen med multipler (near 1,15 til 1,45). Hovedretningen "multipler forklarer mesteparten av near-økningen" ser riktig ut.
- Det som ikke holder: (1) korrelasjonen -0,48/-0,56 mellom underskudd og kopistyrke er delvis mekanisk. I kontrollvinduer under sonen (1580-1880 ms), der det ikke skal være noe underskudd, får jeg -0,25 til -0,38. Med en kopistyrke som ikke bruker de samme near-dataene (nivået på kilden) er sammenhengen -0,16 i sonen og -0,10 til -0,17 i kontrollvinduene, altså ingen forskjell. (2) Test 3 (near mot far) har ingen kraft: nøyaktig samme mønster (near 1,5-1,7 ganger, far 0,6-0,7 ganger mellom dypeste og grunneste fjerdedel, andel 50/50) kommer i alle vinduer, også over og under sonen. (3) "Far er forsinket 6-7 ms bare i sonen" stemmer ikke: forsinkelsen er 4-10 ms hele veien fra 1400 til 1670 ms, med svært lav korrelasjon etter 1490 ms. (4) Å si at dataene "ikke gir støtte for leire med lav Vs" er for sterkt: etter at en tilpasset kopi er trukket fra, står far/near igjen på 0,66-0,72 (kontrollvindu: 0,87-0,89), og i 1470-1510 ms er far lav (0,56-0,6) uten at near er hevet.
- Samtidig kan restunderskuddet ikke tillegges leire. Rett under de sterke reflektorene, i 1354-1366 ms, der verken kopien eller leira (1412-1444 ms fra tid-dybde-kurva) skal ligge, er far/near allerede 0,69-0,75. Og utenfor sonen varierer far/near mellom 0,80 og 1,07 av seg selv. Et rest på 0,7-0,8 er derfor ikke noe diagnostisk.

Min vurdering: dataene peker svakt til moderat mot multipler som hovedforklaring på near-økningen og omtrent halvparten til to tredeler av underskuddet. Resten (far/near 0,7 etter kopifjerning) kan verken tillegges lav Vs eller utelukkes. Test 1 (geometri) har lav, men ikke null kraft hos meg.

## Hovedtall (mine, 94p07 / 94p10; deres i parentes der de finnes)

| Størrelse | 94p07 | 94p10 |
|---|---|---|
| Blokker (200 x 225 m), gyldige | 390 | 364 (deres 382 / 381) |
| far/near sone 1400-1500 / referanse (begge sider) | 0,55 (0,52-0,58) (deres 0,54) | 0,56 (deres 0,54) |
| samme, referanse bare under (1580-1880) | 0,57 (0,54-0,61) (deres 0,57) | 0,59 (deres 0,59) |
| mid/near sone / referanse | 0,69 (0,67-0,71) (deres 0,67) | 0,70 (deres 0,68) |
| kopikoeffisient near / mid / far, forsinkelse | -0,72 / -0,37 / -0,27 ved 110 / 112 / 120 ms | -0,76 / -0,38 / -0,28 ved samme |
| Test 2: r mellom ln(far/near) og ln(kopistyrke) | -0,63 (-0,69 til -0,56); referanse bare under -0,52 | -0,66; bare under -0,52 |
| Test 2: samme målt med bare formen på kopien | -0,40 | -0,50 |
| Test 2: kontrollvinduer 1580-1880, samme test | -0,25 til -0,35 | -0,24 til -0,40 |
| Test 2: uavhengig prediktor (nivå på kilden), sone / kontroll | -0,16 / -0,10 til -0,16 | -0,16 / -0,11 til -0,17 |
| far/near sone/ref, svakeste tredel av kopien (uavhengig halvdel) | 0,69-0,70 | 0,73 |
| samme, svakeste tidel (uavhengig halvdel) | 0,77-0,85 | 0,80-0,83 |
| Nivå i sonen mot referanse under, hele området: near / mid / far | 1,45 / 1,0 / 0,8-0,9 (30 ms vinduer 1414-1444 ms: near 1,7-2,0, far 0,9-1,1) | tilsvarende |
| far/near etter fjerning av tilpasset kopi (110 og 220 ms), sone | 0,67-0,72 (rå 0,57) | 0,68-0,74 (rå 0,59) |
| samme i kontrollvindu 1580-1680 | 0,87-0,94 (rå 0,86) | 0,87-0,94 (rå 0,86) |
| Test 1: sd av avstanden 1316 til 1386 ms | 4,2 ms (5-95 %: 12,3 ms); i 300 m-blokker 3,2 ms | 4,6 ms (15,0 ms) |
| Test 3: dypeste/grunneste fjerdedel, near / mid / far (same data) | 1,54 / 1,19 / 0,61 | 1,47 / 1,16 / 0,59 |
| Test 3: samme i vinduet 1580-1680 (kontroll) | 1,60 / 1,05 / 0,60 | 1,62 / 1,06 / 0,60 |
| Andel av variasjonen fra far / near | 0,53 / 0,47 i sonen, 0,47-0,54 i andre vinduer | 0,56 / 0,44 |

## Hva som er gjort

1. Egen blokkoppbygging direkte fra kubene (mmap, én stripe om gangen, under 800 MB): celler på 8 x 6 spor (100 x 75 m), tid 1000-1950 ms, summert til blokker på 200 x 225 m (og 300 m-blokker som kontroll). Levende spor og dekningen for 94p10 er håndtert (98 % dekning). Tidsaksen er sjekket mot havbunnen: sterkeste utslag på near/ful ved 106-108 ms, i samsvar med 81 m vanndybde. Brønnen 15/9-13 sitt tid-dybde-forhold gir leira (1383-1415 m MD) ved 1412-1444 ms.
2. Hver blokktrase er flatet ved krysskorrelasjon av hele stakken i 1236-1346 ms mot brønnblokka. Energi er målt som RMS av blokkens midlere trase og som produkt av to rutemønsterhalvdeler (for å fjerne støy). Dette er ikke de plukkede hendelsene den andre analysen brukte.
3. Kopi: minste kvadraters koeffisient mellom near(t) og near(t - lag) for pakken 1256-1356 ms, lag søkt 90-134 ms, per stakk.
4. Usikkerhet: bootstrap over felt på 3 x 4 blokker (600 x 900 m), og korrelasjon på gitter med hver tredje blokk. Rangkorrelasjon er rapportert ved siden av.
5. Kontroller jeg la til: (a) kontrollvinduer 1580-1880 ms med samme test, (b) prediktor for kopistyrke som ikke bruker sonens egne near-data, (c) kopifjerning i vindu og referanse, (d) fin tidsprofil i 30 ms-vinduer, (e) forsinkelse av far mot near langs hele tidsaksen.

## Hva som er usikkert eller svakt

- Test 2 deler near-data mellom underskuddet og kopistyrken. Kontrollvinduene viser at sammenhengen er -0,25 til -0,38 uten noen forventet effekt, så bare overskuddet (omtrent -0,25 til -0,3) er tegn på noe ekte. Den uavhengige prediktoren gir ingen forskjell mellom sone og kontroll. Det er likevel ikke et bevis mot multipler: nivået på kilden er et svakt mål på kopistyrke (r med near i sonen bare 0,25), og kilden i kontrollvinduene er svakere.
- Min test 1 er ikke et fullverdig svar: jeg har bare standardavviket til avstanden (4,2 ms, spenn 12 ms). Jeg har ikke gjort stigningstallstesten (følger minimum de sterke reflektorene eller sonen). Notatets 2,5 ms og mitt 4,2 ms skyldes ulik flating, og spennet over 8 ms ligger akkurat over terskelen. Testen er svak, ikke uten kraft.
- Kopifjerning med minste kvadrater fjerner også ekte primærenergi som tilfeldigvis ligner på signalet 110 ms over, så 0,67-0,72 er et nedre anslag for hva multipler forklarer. Det er ikke en ren måling av restunderskuddet.
- Underskuddet starter allerede ved 1348-1366 ms (far/near 0,69-0,75, near 0,85-0,87, far 0,6-0,64), over både kopiene og leira. Det tyder på at bølgeform og avstandsavhengig haleeffekt under en sterk reflektorpakke gir en del av effekten uavhengig av begge forklaringene.
- Utenfor 1380-1560 ms varierer far/near mellom 0,80 og 1,07. Et restunderskudd på 20-30 % er derfor i samme størrelse som vanlig variasjon.
- Sonen ble valgt fordi den er lav (utvalgseffekt). Kontrollvinduene er likevel valgt uten å se på resultatene.
- 94p10 er samme feltdata som 94p07 og er ikke en uavhengig bekreftelse (underskuddet korrelerer 0,91 mellom dem hos dem).
- Jeg har ikke etterprøvd modellen eller prediksjoner.json.

## Hva som ville endret min vurdering

Gathers som viser om hendelsene i sonen flater ut som primærer, multippelfjerning på gathers, eller en skjærlogg gjennom laget.
