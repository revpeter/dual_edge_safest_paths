# Results description
**Summary táblák oszlopai**
- network - Hálózat neve
- fugg_jobb/dual_jobb - A duális szemüvegen keresztül az adott módszer hányszor talált jobb útvonalat.
- eq - A két módszer ennyiszer ugyanolyan útvonalat talált a duális szemüvegen keresztül nézve.
- eq/fugg/dual_meanLen - Átlagos hossza az útvonalaklak.
- fuggBet/dualBet_meanEdAvb - Azoknak az útvonalaknak az átlagos él-duális féle availability-je amiknél az adott módszer jobb utat adott.
- fugg/dual_meanZpct - Az útvonalakban a 0-as edge-k átlagos aránya az útvona hosszához.
- th - Az srlg sorbarendezéssel talált threshold az adott hálózatban. Ha a th_mod értéke 1 akkor volt vágás a thresholdnál, ha 0 akkor nem volt (a köztes értékek eredményeit nem töltöttem fel).

**simResult táblák oszlopai**
- simID - Az adott szimuláció id-ja.
- start - A node ahonnan indulunk.
- target - A node ahova el akarunk jutni.
- len_sht/fugg/djk - A legrövidebb/független/él-duális módszer által talált út hossza.
- avb_sht/fugg/djk - A független módszer szerinti  availability a legrövidebb/független/él-duális útvonalakra.
- ed_avb_sht/fugg/djk - Az él-duális féle availability a legrövidebb/független/él-duális útvonalakra.
- ed_avb_diff - ed_avb_djk - ed_avb_fugg. Ha >0 akkor a független módszer útvonala a jobb. Ha <0 akkor az él-duális módszer útvonala a jobb.
- pct_zero_djk/fugg/sht - Az adott útvonalban a 0-s edge-k arány az útvonal hosszához.
- 
