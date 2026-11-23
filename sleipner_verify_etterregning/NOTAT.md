# Etterregning av Sleipner-testen (uavhengig kontroll)

Kjørt 19. november 2026, ca. 22.05–22.45, på den leide maskinen og Linux-PC-en. Egen behandling og egne mål fra rå-gatherne (`gather_G_*`, `gather_FS_*`). `slp.py` og `s3_stack_measure.py` er ikke importert. Alle tall under kommer fra skriptene i denne mappen.

## Dom

Tallene i rapporten lar seg regne etter. Konklusjon 2 og 3 holder. Konklusjon 1 holder som tall, men «gir mønsteret» er for sterkt. Konklusjon 4 holder og bør stå tydeligere.

## Hva jeg fikk

| Mål (sone/referanse for far/near, 25 Hz, 48° mute) | Uten havflate: trend / −7 / −15 / −20 % | Med havflate: trend / −7 / −15 / −20 % |
|---|---|---|
| RMS i vinduer | 1,00 / 0,97 / 0,98 / 1,00 | 0,62 / 0,57 / 0,53 / 0,53 |
| Topp per reflektor (3 i sonen, 5 utenfor) | 0,91 / 1,00 / 1,07 / 1,04 | 0,47 / 0,43 / 0,40 / 0,38 |
| Omhylningstopp per reflektor | 0,82 / 0,77 / 0,72 / 0,82 | 0,53 / 0,51 / 0,45 / 0,48 |

- Rapporten har RMS 1,03 uten og 0,52 med havflate. Mine stakker korrelerer 0,9995 (near), 0,9999 (mid) og 0,97–0,98 (far) med rapportens.
- mid/near faller også med havflate: 0,43 (RMS), 0,27–0,61 (per reflektor).
- Med 50 % strekk-mute: 1,05 uten, 0,48 med. Med 40° mute: 0,81 uten, 0,50 med (referansen krymper til 459 av 620 ms).

## Det som bør rettes i rapporten

1. **Underskuddet er ikke spesielt for sonen.** Samme mål i alle 84 ms-vinduer mellom 1,1 og 1,9 s: uten havflate 0,52–1,59 (10–90 %: 0,64–1,34), og 12 % av vinduene ligger under 0,65. Med havflate 0,37–1,49, og 26 % ligger under 0,65. Sonen (0,62) er på 19-prosentilen. Per reflektor spenner far/near 0,38–1,89 utenfor sonen uten havflate. I de ekte dataene er spennet 0,38–0,68.
2. **I modellen er det mid som faller mest, og near som øker.** Total/primær (RMS) i sonen mot utenfor: near 1,83 mot 1,46, mid 0,71 mot 1,26, far 0,94 mot 1,21. De ekte dataene har normal mid/near med hovedmetoden (0,94–1,00).
3. **Forsinkelsen på kopiene i de ekte dataene passer ikke med offsetområdene.** Etter NMO er forsinkelsen 109,6 ms ved 150 m, 111,9 ved 760 m, 122,3 ved 1760 m og 145 ved 3000 m. Målt 110 / 112 / 119 ms svarer til ca. 360 / 780 / 1500 m, og andre kopi (220 / 224 / 234 ms) til ca. 370 / 800 / 1380 m. En far-stakk fra 1760 m og utover kan ikke gi 119 ms.
4. **Alle referansereflektorene ligger over sonen** (1109–1331 ms hos meg, 1121–1347 ms i rapporten). Ingen ligger under 1560 ms.
5. **Muten i rapporten er hakkete.** Antall far-spor i sonen hopper mellom 24 og 51 fra prøve til prøve (159 fall mellom 1,1 og 1,9 s). Med jevn mute får jeg 0,62 i stedet for 0,52 med havflate. Konklusjonen endres ikke.
6. **Near-kopien er målt for sterk.** Jeg måler −0,91 og −1,00 ved 111 ms (rapporten −0,93). Teori for én reflektor er −0,75 ved 111 ms. På valideringsmodellen måler jeg −0,74 til −0,76 mot teori −0,757 til −0,752.
7. **Linjen om 1 m rutenett er fortsatt en plassholder.** 1 m-kjøringen ble ferdig ca. 22.31. Med min behandling (`e15`): stakkene korrelerer 0,997 (near), 0,997 (mid) og 0,996 (far) med 2 m, og sone/referanse med havflate er 0,61 mot 0,62 (RMS), 0,47 mot 0,47 (topp) og 0,53 mot 0,53 (omhylning). Amplituden i 1 m-kjøringen er 0,249 av 2 m-kjøringen (kilden skalerer med cellearealet), så rå amplituder kan ikke sammenlignes direkte.
8. **Trykket er ikke null i havflaten i de elastiske speilkjøringene.** Rapporten sier at det er det. Se neste avsnitt. Virkningen på multiplene er liten.
9. **«Aldri underskudd uten havflate» er for sterkt.** Med 40° mute får jeg 0,80–0,82 uten havflate, og omhylningsmålet gir 0,72–0,82 med 48°.

## Havflaten: speilet er ikke helt symmetrisk

- Jeg kjørte rapportens egen `s2_fd.py` på valideringsmodellen med ekstra mottakere i z = 0 (`e6`). Trykket i z = 0 er opptil 24 % av trykket på 8 m på det verste sporet (1,1 % av største utslag i hele opptaket).
- Årsak (`e10`–`e12`): Devito 4.8.23 midler skjærmodulen til de forskjøvne t_xz-punktene med `safe_harmonic`, og vernet mot null ser bare på noden (x, y). Ved den ekte havbunnen blir skjærmodulen 0 (glidning, riktig). Ved speilbildet av havbunnen blir den 2 × sedimentets. Speilmodellen er dermed ikke et eksakt speilbilde.
- Med bare vann, eller med væske under havbunnen, er symmetrien eksakt (1e−7 til 1e−5). Med elastisk havbunn er trykket i z = 0 7–16 % (RMS) av trykket på 8 m.
- Hvor mye betyr det (`e13`, `e14`, liten modell, reflektor på 400 m, vinkel i vannet 0–34°)? Jeg laget en variant som er eksakt symmetrisk og har glidning på begge havbunner (trykk i z = 0: 1e−6). Første kopi/primær er lik innen 0,004 i de to (0,988 mot 0,987 ved 0°, 0,813 mot 0,811 ved 34°). Multiplene avviker 1–14 % (RMS av forskjellen), nivået innen 2 %, korrelasjon ≥ 0,990. Primæren med spøkelser avviker 0,05 % ved 0° og 3,3 % ved 34°.
- Konklusjonene berøres ikke, men setningen «Da er trykket null i havflaten» bør rettes.

## Kontroll av delene

| Del | Resultat |
|---|---|
| Havflaten, bare vann | Vann-kjøringen mot eksakt løsning (linjekilde med speilkilde): amplitude innenfor −12 til +5 % ut til 3150 m, korrelasjon 0,968–0,9999. En refleksjonskoeffisient på −0,99 ville gitt avvik på ca. 100–230 % fra 1500 m. |
| Havflaten, elastisk havbunn | Ikke eksakt null, se avsnittet over. Første kopi på valideringsmodellen stemmer likevel: −0,739 til −0,759 mot teori −0,757 til −0,752, forsinkelse innen 0,5 ms. |
| Primær lik med og uten havflate | Valideringsmodellen: amplitudeforhold 0,9997 ved 0 m, 0,9918 ved 800 m. |
| Tid null | Avvik −0,05 ms ved 150 m mot eksakt løsning. Havbunnsrefleksen ved 150 m: 149,0 ms, ventet 149,2. |
| 35 → 25 Hz | Eksakt i form: amplitudeforhold 0,9997–1,03, korrelasjon ≥ 0,991. Rapportens filter har en konstant faktor 25/35 som ikke påvirker forhold. |
| NMO | Rest høyst 1,3 ms for primærene ut til 2250 m. |
| Klipping | Alle gathere endelige. Ingen avtrapping. |
| Egen Zoeppritz | Energibevaring 1,3e−15 (fast–fast) og 1,8e−15 (væske–fast), 2000 tilfeldige tilfeller hver. |

## Vs i laget, bare primærer (egen Zoeppritz)

- Bunnen av laget alene: far/near 1,12 → 0,91 → 0,69 → 0,57 for 0, 7, 15, 20 % lavere Vs (rapporten: 1,10 → 0,68 ved 15 %).
- Toppen får samtidig sterkere far: Rpp på far +0,033 ved toppen og −0,028 ved bunnen for 15 %.
- Som bygget, per reflektor ved 15 %: 1410 ms 1,22 → 1,03, 1434 ms 1,28 → 1,17, 1458 ms 1,32 → 1,61. Sone/referanse (RMS): 0,975 → 0,969 → 0,986 → 1,006.
- FD mot Zoeppritz for Vs-virkningen: styrke 1,02–1,03 (near), 1,00 (mid), 1,16–1,18 (far).

## Trepunktsfilter (eget)

Etter filter: 0,74–1,36 (trend) og 0,90–1,68 (−15 %), der fasiten er 0,72–1,07. Filteret bommer oftest oppover. Ingen av mine filtrerte verdier er så lave som 0,62–0,65, som står igjen i de ekte dataene. På far fjerner det ingenting (rest 0,67 mot 0,62 før).

## Ikke kontrollert

Variantene med havbunn-R 0,37 og 0,43, Vs −30 og −40 %, bindingen til den ekte trasen, og selve målingene på de ekte dataene.

## Filer

| Fil | Innhold |
|---|---|
| `e1_lib.py` | Egen behandling og egne mål |
| `e2_main.py`, `e2.log`, `e2_results.json`, `e2_stacks.npz` | Stakker, sone/referanse, kopier, NMO, fold |
| `e3_zoep.py` | Egen Zoeppritz med energitest |
| `e4_surface_analytic.py`, `e4.log`, `e4_results.json` | Vann-kjøringene mot eksakt løsning, 35 → 25 Hz |
| `e5_vs_zoep.py`, `e5.log`, `e5_results.json`, `e5_zoep_stacks.npz` | Vs i laget, bare primærer |
| `e6_make.py`, `e6_fs_surface.py`, `e6.log`, `e6run/` | Rapportens `s2_fd.py` med ekstra mottakere i z = 0 |
| `e7_extras.py`, `e7.log`, `e7_results.json` | Alle vinduer, valideringskopier, forsinkelse mot offset, sammenligning med rapportens stakker, filter |
| `e8_crosscheck.py`, `e8_results.json` | Væske–fast energitest, FD mot Zoeppritz |
| `e9_primary_fs_vs_g.py`, `e9_results.json` | Primær med og uten havflate på alle offset |
| `e10_mirror_symmetry.py`, `e11_stencils.py`, `e12_where_asym.py` med logg og json | Hvorfor speilet ikke er symmetrisk |
| `e13_impact.py`, `e14_symslip.py` med logg, json og npz | Hvor mye det betyr for multiplene |
| `e15_dx1.py`, `e15.log`, `e15_results.json` | 1 m mot 2 m rutenett |
