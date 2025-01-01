import os
import matplotlib.pyplot as plt
from typing import List, Tuple
import math
import datetime
from dataclasses import dataclass
#from jinja2 import Environment, FileSystemLoader
import bisect
#from weasyprint import HTML
#import pdfkit

base_path_reiberschema = os.path.abspath(os.path.join(os.path.dirname(__file__)))


"""
Quellen für die Formeln:

https://www.horeiber.de/?page_id=306

https://www.horeiber.de/pdf/eins.pdf

Die intrathekale Fraktion, IgIF(%) ist direkt von den
Prozentlinien der Quotientendiagramme ablesbar (Abb.
46-2). Die Prozentlinie der Quotientendiagramme in
Abb. 46-2, z. B. für 20%, wird nach der Formel mit
(1 – QLim/QIgG) = 0,2 oder QIgG = 1,25 × QLim be-
rechnet  (eins S.36)


Reiber H (1994). The hyperbolic function: a mathematical solution of the protein flux/CSF flow model 
for blood-CSF barrier function J Neurol Sci 126:243-245.

Reiber H (1994). Flow rate of cerebrospinal fluid (CSF)- a concept common to normal blood-CSF barrier function 
and to dysfunction in neurological diseases. J Neurol Sci 122:189-203
"""

"""
Testdaten aus csv Datei:

20438528|Petzolt|Bernd|08.05.1967|M|C110|Kreisklinikum Siegen GmbH|Station 31 Neurologie|02717051850|57076|SIEGEN|Weidenauerstr. 76||6920438528|21.10.2023|61|6.61|mg/dl|67|49.10|mg/dl|
20438528|Petzolt|Bernd|08.05.1967|M|C110|Kreisklinikum Siegen GmbH|Station 31 Neurologie|02717051850|57076|SIEGEN|Weidenauerstr. 76||5120438528|21.10.2023|1|1160.00|mg/dl|7|4030.00|mg/dl|

tnr|Nachname|Vorname|Gebdat|Geschlecht|?|Einsender|Station|Tel|PLZ|Ort|Adresse|-|MatNummer|Liste: Code|Ergebnis|Einheit

Verfahrenscodes:
1 => IGGLS
7 => BALBN 
61 => IGGL
67 => ALBL
"""

auftrag_map = dict()


@dataclass(frozen=True)
class Einsender:
    name: str
    station: str
    adresse: str
    ort: str
    telefon: str


@dataclass(frozen=True)
class Patient:
    vorname: str
    nachname: str
    geburtsdatum: str  # Format %d.%m.%Y
    geschlecht: str

    @property
    def pat_anzeige(self):
        if self.geschlecht == 'F':
            return 'Patientin'
        return 'Patient'

    @property
    def alter(self) -> int:
        now = datetime.datetime.now()
        geb_dat_dt = datetime.datetime.strptime(self.geburtsdatum, '%d.%m.%Y')
        delta = now - geb_dat_dt
        return int(delta.days / 365.25)


@dataclass(frozen=True)
class AnalytIdentifikation:
    verf_liquor: str
    verf_serum: str
    code_liquor: str
    code_serum: str


class AnalytErgebnis:

    def __init__(self, bezeichnung: str, identifikation: AnalytIdentifikation):
        self.bezeichnung = bezeichnung
        self.identifikation = identifikation
        self.wert_liquor = None
        self.liquor_einheit = ''
        self.wert_serum = None
        self.serum_einheit = ''
        self.quotient: float = 0.0
        self.synthese_pct = 0.0  # Fraktion intrathekale Synthese
        self.bild_datei = ''

    @property
    def quotient_string(self):
        return f'{self.quotient:.2f}'

    @property
    def liquor_ergebnis_tabelle(self) -> str:
        if not self.ergebnis_ist_vorhanden():
            return ''
        return f'{float(self.wert_liquor):.1f} {self.liquor_einheit}'

    @property
    def serum_ergebnis_tabelle(self) -> str:
        if not self.ergebnis_ist_vorhanden():
            return ''
        return f'{float(self.wert_serum):.1f} {self.serum_einheit}'

    @property
    def quotient_ergebnis_tabelle(self):
        if not self.ergebnis_ist_vorhanden():
            return ''
        quotient = self.quotient * 1000
        return f'{quotient:.1f}'

    @property
    def synthese_anzeige(self):
        if not self.ergebnis_ist_vorhanden():
            return ''
        if self.synthese_pct == 0:
            return '0%'
        return f'{self.synthese_pct*100:.1f}%'

    def ergebnis_ist_vorhanden(self) -> bool:
        return self.wert_liquor is not None

    def berechne_quotient(self):
        if self.wert_liquor and self.wert_serum:
            self.quotient = float(self.wert_liquor) / float(self.wert_serum)


class Ergebnisse:

    def __init__(self, liste_aus_csv: List[str]):
        self.albumin = AnalytErgebnis('Albumin', AnalytIdentifikation('ALBL', 'BALBN', '67', '7'))
        self.igg = AnalytErgebnis('IgG', AnalytIdentifikation('IGGL', 'IGGLS', '61', '1'))
        self.iga = AnalytErgebnis('IgA', AnalytIdentifikation('IGAL', 'IGALS', '62', '2'))
        self.igm = AnalytErgebnis('IgM', AnalytIdentifikation('IGML', 'IGMLS', '63', '3'))
        self.werte_liste_aus(liste_aus_csv)

    def werte_liste_aus(self, liste_aus_csv: List[str]):
        def setze_wert_und_einheiten(erg: AnalytErgebnis, liquor_identifier_liste, serum_identifier_liste):
            if bezeichnung in liquor_identifier_liste:
                erg.wert_liquor = wert
                erg.liquor_einheit = einheit
            elif bezeichnung in serum_identifier_liste:
                erg.wert_serum = wert
                erg.serum_einheit = einheit

        ergebnisse = [liste_aus_csv[x:x + 100] for x in range(0, len(liste_aus_csv), 3)]
        for ergebnis in ergebnisse:
            bezeichnung = ergebnis[0]
            wert = ergebnis[1]
            einheit = ergebnis[2]

            setze_wert_und_einheiten(self.albumin, [self.albumin.identifikation.verf_liquor, self.albumin.identifikation.code_liquor],
                                     [self.albumin.identifikation.verf_serum, self.albumin.identifikation.code_serum])
            setze_wert_und_einheiten(self.igg, [self.igg.identifikation.verf_liquor, self.igg.identifikation.code_liquor],
                                     [self.igg.identifikation.verf_serum, self.igg.identifikation.code_serum])
            setze_wert_und_einheiten(self.iga, [self.iga.identifikation.verf_liquor, self.iga.identifikation.code_liquor],
                                     [self.iga.identifikation.verf_serum, self.iga.identifikation.code_serum])
            setze_wert_und_einheiten(self.igm, [self.igm.identifikation.verf_liquor, self.igm.identifikation.code_liquor],
                                     [self.igm.identifikation.verf_serum, self.igm.identifikation.code_serum])

        self.albumin.berechne_quotient()
        self.igg.berechne_quotient()
        self.iga.berechne_quotient()
        self.igm.berechne_quotient()

    def als_liste(self) -> List[AnalytErgebnis]:
        return [self.igg, self.iga, self.igm]


class Probe:

    def __init__(self, proben_nr, material, datum):
        self.proben_nr = proben_nr
        self.material = material
        self.datum = datum


class CsvZeile:

    """
    #tnr|Nachname|Vorname|Gebdat|Geschlecht|?|Einsender|Station|Tel|PLZ|Ort|Adresse|-|MatNummer|
    Liste: Code|Ergebnis|Einheit
    """

    def __init__(self, csv_zeile: str):
        csv_zeile = csv_zeile.strip()
        if csv_zeile[-1] == '|':
            csv_zeile = csv_zeile[:-1]
        data = csv_zeile.split('|')
        self.auftrag_nr = data[0]
        self.pat_nachname = data[1]
        self.pat_vorname = data[2]
        self.pat_geburtsdatum = data[3]
        self.pat_geschlecht = data[4]
        self.eins_name = data[6]
        self.eins_station = data[7]
        self.eins_tel = data[8]
        self.eins_plz = data[9]
        self.eins_ort = data[10]
        self.eins_adre = data[11]
        self.barcode = data[13]
        self.proben_datum = data[14]
        self.ergebnisse = data[15:]


class Auftrag:

    def __init__(self, csv_zeile: str):
        data = CsvZeile(csv_zeile)
        self.datum_heute = datetime.datetime.now().strftime('%d.%m.%Y')
        self.auftrag_nr = data.auftrag_nr
        self.patient = Patient(data.pat_vorname, data.pat_nachname, data.pat_geburtsdatum, data.pat_geschlecht)
        ort = f'{data.eins_plz} {data.eins_ort}'
        self.einsender = Einsender(data.eins_name, data.eins_station, data.eins_adre, ort, data.eins_tel)
        verfahren_liste = data.ergebnisse
        barcode = data.barcode
        self.proben: List[Probe] = []  # Barcodes mit Materialindex
        self.neue_probe(barcode, data.proben_datum)
        self.ergebnisse: Ergebnisse = Ergebnisse(verfahren_liste)
        self.delpq: str = '0.0'  # Verfahren DELPQ

    def werte_csv_zeile_aus(self, csv_zeile: str):
        data = CsvZeile(csv_zeile)
        self.neue_probe(data.barcode, data.proben_datum)
        verfahren_liste = data.ergebnisse
        self.ergebnisse.werte_liste_aus(verfahren_liste)

    def neue_probe(self, barcode: str, datum: str):
        if barcode.startswith('69'):
            self.proben.append(Probe(barcode, 'Liquor', datum))
        elif barcode.startswith('51'):
            self.proben.append(Probe(barcode, 'Serum', datum))

    def berechne_delpq(self):
        """
        Delpech-Lichtblauquotient = Qigg/Qalb
        Ist der Wert des Delpech-Lichtblau-Quotienten > 0,7 spricht
        dies für intrathekale IgG-Synthese
        """
        delpq = self.ergebnisse.igg.quotient / self.ergebnisse.albumin.quotient
        print("delpq", delpq)
        self.delpq = f'{delpq:.2f}'

    @property
    def albumin_ergebnis(self) -> AnalytErgebnis:
        return self.ergebnisse.albumin

    @property
    def igg_datei(self) -> str:
        return self.ergebnisse.igg.bild_datei

    @property
    def iga_datei(self) -> str:
        return self.ergebnisse.iga.bild_datei

    @property
    def igm_datei(self) -> str:
        return self.ergebnisse.igm.bild_datei


def _create_orders_from_csv(file_path: str):
    with open(file_path, 'r') as csv_datei:
        lines = csv_datei.readlines()
        for line in lines:
            data = line.split('|')
            auftrag_nr = data[0]
            if auftrag_nr in auftrag_map:
                auftrag_map[auftrag_nr].werte_csv_zeile_aus(line)
            else:
                auftrag_map[auftrag_nr] = Auftrag(line)

    for auftrag_nr, auftrag in auftrag_map.items():
        auftrag.berechne_delpq()


class DiagramDimension:
    x_values = [x * 0.0001 for x in range(10, 1000)]
    y_max = 0.15


class CurveParameters:

    """
    Parameter fuer ein bestimmtes Immunoglobulin basierend auf:
    Reiber H (1994). Flow rate of cerebrospinal fluid (CSF)- a concept common to normal blood-CSF barrier function
    and to dysfunction in neurological diseases. J Neurol Sci 122:189-203
    """

    def __init__(self, name: str, ab_lower: float, b2_lower: float, c_lower: float, ab_upper: float, b2_upper: float, c_upper: float):
        self.name = name
        self.ab_lower = ab_lower
        self.b2_lower = b2_lower * 0.000001
        self.c_lower = c_lower * 0.001
        self.ab_upper = ab_upper
        self.b2_upper = b2_upper * 0.000001
        self.c_upper = c_upper * 0.001


class ImmunglobulineCurves:

    """
    Haelt die Listen aus x und y Werten die im Diagramm fuer die
    Kurven verwendet werden
    """

    def __init__(self, params: CurveParameters):
        self.params = params
        self.x = DiagramDimension.x_values
        self.upper: List[float] = []
        self.lower: List[float] = []
        self.pct_20: List[float] = []
        self.pct_40: List[float] = []
        self.pct_60: List[float] = []
        self.pct_80: List[float] = []
        self.calculate_data_points()

    @property
    def name(self):
        return self.params.name

    def last_visible_20pct_x(self, y_max) -> Tuple[float, float]:
        next_index = bisect.bisect_left(self.pct_20, y_max)
        next_index = next_index - 1
        return self.x[next_index], self.pct_20[next_index]

    def last_visible_40pct_x(self, y_max) -> Tuple[float, float]:
        next_index = bisect.bisect_left(self.pct_40, y_max)
        next_index = next_index - 1
        return self.x[next_index], self.pct_40[next_index]

    def last_visible_60pct_x(self, y_max) -> Tuple[float, float]:
        next_index = bisect.bisect_left(self.pct_60, y_max)
        next_index = next_index - 1
        return self.x[next_index], self.pct_60[next_index]

    def last_visible_80pct_x(self, y_max) -> Tuple[float, float]:
        next_index = bisect.bisect_left(self.pct_80, y_max)
        next_index = next_index - 1
        return self.x[next_index], self.pct_80[next_index]

    def calculate_data_points(self):
        for alb_q in self.x:
            lower_y = self.get_lower_lim_at_coordinate(alb_q)
            self.lower.append(lower_y)
            upper_y = self.get_upper_lim_at_coordinate(alb_q)
            self.upper.append(upper_y)

            pct_20 = self.get_upper_lim_at_coordinate(alb_q) * 1.25
            self.pct_20.append(pct_20)
            pct_40 = self.get_upper_lim_at_coordinate(alb_q) * (10/6)
            self.pct_40.append(pct_40)
            pct_60 = self.get_upper_lim_at_coordinate(alb_q) * 2.5
            self.pct_60.append(pct_60)
            pct_80 = self.get_upper_lim_at_coordinate(alb_q) * 5
            self.pct_80.append(pct_80)

    def get_lower_lim_at_coordinate(self, q_albumin: float) -> float:
        return self.params.ab_lower * math.sqrt(pow(q_albumin, 2) + self.params.b2_lower) - self.params.c_lower

    def get_upper_lim_at_coordinate(self, q_albumin: float, factor=1.0) -> float:
        return self.params.ab_upper * math.sqrt(pow(q_albumin*factor, 2) + self.params.b2_upper) - self.params.c_upper

    def setze_synthese_pct_fuer_ergebnis(self, ergebnis: AnalytErgebnis, albumin_ergebnis: AnalytErgebnis):
        """
        pct = (Q_igx-Q_lim) / Q_igx
        :param ergebnis: Ergebnis der IgX Messung
        :param albumin_ergebnis: Ergebnis der Albumin Messung fuer die Berechnung von Q_lim
        """
        if not ergebnis.ergebnis_ist_vorhanden():
            return
        q_lim = self.get_upper_lim_at_coordinate(albumin_ergebnis.quotient)
        syn_pct = (ergebnis.quotient - q_lim) / ergebnis.quotient
        if syn_pct < 0:
            syn_pct = 0
        ergebnis.synthese_pct = syn_pct


def plot_aufrag_images(auftrag: Auftrag):
    immun_globulin_param_map = dict(
    IgG=CurveParameters('IgG', 0.33, 2, 0.3, 0.93, 6, 1.7),
    IgA=CurveParameters('IgA', 0.17, 74, 1.3, 0.77, 23, 3.1),
    IgM=CurveParameters('IgM', 0.04, 442, 0.82, 0.67, 120, 7.1)
    )

    age = auftrag.patient.alter  # QAlb = (4 + Alter(Jahre)/15) × 10–3
    print("age", age)
    q_alb = (4 + age/15) * 0.001

    y_max = 0.15

    # Schleife ueber alle Immunglobuline:
    for analyt_ergebnis in auftrag.ergebnisse.als_liste():

        ig_select = analyt_ergebnis.bezeichnung

        if not analyt_ergebnis.ergebnis_ist_vorhanden():
            datei_name = f'leer/{ig_select}_leer.svg'
            print("Datei", datei_name)
            analyt_ergebnis.bild_datei = datei_name
            continue

        igs_curves = ImmunglobulineCurves(immun_globulin_param_map[ig_select])

        igs_curves.setze_synthese_pct_fuer_ergebnis(analyt_ergebnis, auftrag.albumin_ergebnis)
        print("synthese", analyt_ergebnis.synthese_pct)

        fig = plt.figure()

        ax = fig.add_subplot(1, 1, 1)
        ax.plot(igs_curves.x, igs_curves.lower, color=(0, 0, 0))
        ax.set_xscale('log')
        ax.set_yscale('log')
        # Referenzbereich:
        ax.plot(igs_curves.x, igs_curves.lower, color=(0, 0, 0), linewidth=3)
        ax.plot(igs_curves.x, igs_curves.upper, color=(0, 0, 1), linewidth=3)
        # Prozentlinien intrathekale Fraktion:
        last_20pct_datapoint = igs_curves.last_visible_20pct_x(y_max)
        last_40pct_datapoint = igs_curves.last_visible_40pct_x(last_20pct_datapoint[1])
        last_60pct_datapoint = igs_curves.last_visible_60pct_x(last_20pct_datapoint[1])
        last_80pct_datapoint = igs_curves.last_visible_80pct_x(last_20pct_datapoint[1])

        y_pct_lines = last_20pct_datapoint[1]

        ax.plot(igs_curves.x, igs_curves.pct_20, color=(0, 0, 0, 0.5), linestyle='dashed')
        ax.annotate(xy=(last_20pct_datapoint[0], y_pct_lines), text='20%', weight='bold')

        ax.plot(igs_curves.x, igs_curves.pct_40, color=(0, 0, 0, 0.5), linestyle='dashed')
        ax.annotate(xy=(last_40pct_datapoint[0], y_pct_lines), text='40%', weight='bold')

        ax.plot(igs_curves.x, igs_curves.pct_60, color=(0, 0, 0, 0.5), linestyle='dashed')
        ax.annotate(xy=(last_60pct_datapoint[0], y_pct_lines), text='60%', weight='bold')

        ax.plot(igs_curves.x, igs_curves.pct_80, color=(0, 0, 0, 0.5), linestyle='dashed')
        ax.annotate(xy=(last_80pct_datapoint[0], y_pct_lines), text='80%', weight='bold')

        ax.grid()
        grid_lines = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]
        ax.set_xticks(grid_lines)
        ax.set_yticks(grid_lines)
        x_labels = [1, 2, 5, 10, 20, 50, '100']
        ax.set_xticklabels(x_labels, weight='bold')
        y_labels = [1, 2, 5, 10, 20, 50, '100']
        ax.set_yticklabels(y_labels, weight='bold')

        # Messpunkt markieren:
        if analyt_ergebnis.ergebnis_ist_vorhanden():
            plt.scatter(auftrag.ergebnisse.albumin.quotient, auftrag.ergebnisse.igg.quotient, c='black', marker='D', s=100)

        # Axen:
        ax.set_xlim([0.001, 0.1])
        ax.set_ylim([0, y_max])

        ax.tick_params(axis="x", direction="in", pad=-20)
        ax.tick_params(axis="y", direction="in", pad=-30)

        ax.text(0.01, 0.75, '$Q_{' + igs_curves.name + '}$', size=22, transform=ax.transAxes, fontweight='bold')
        ax.text(0.15, 0.75, '(x $10^{-3}$)', size=12, transform=ax.transAxes)

        ax.text(0.7, 0.1, '$Q_{Alb}$', size=22, transform=ax.transAxes, fontweight='bold')
        ax.text(0.85, 0.1, '(x $10^{-3}$)', size=12, transform=ax.transAxes)

        plt.vlines(x=q_alb, ymin=igs_curves.get_lower_lim_at_coordinate(q_alb), ymax=y_max,
                   colors='purple',
                   label='vline_multiple - full height')

        #plt.show()
        datei_name = f'tmp/{auftrag.auftrag_nr}_{ig_select}.svg'
        analyt_ergebnis.bild_datei = datei_name
        print("Datei", datei_name)
        plt.savefig(datei_name, format='svg',  bbox_inches='tight')


def create_images_for_file(file_path: str):
    _create_orders_from_csv()
    for auftrag_nr, auftrag in auftrag_map.items():
        plot_aufrag_images(auftrag)

