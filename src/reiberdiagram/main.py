import os
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional
import math
import datetime
from dataclasses import dataclass
import bisect

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
auftrag_list = []


class Immunglobulin:
    IGA = "IgA"
    IGG = "IgG"
    IGM = "IgM"


class ImageType:
    PNG = "png"
    JPG = "jpg"
    SVG = "svg"


class AnalytResult:

    def __init__(self, name: Immunglobulin, value_serum: float, value_csf: float):
        self.name = name
        self.value_csf = value_csf
        self.value_serum = value_serum
        self.quotient: float = 0.0
        self.synthese_pct = 0.0  # Fraktion intrathekale Synthese
        self.bild_datei = ''

    @property
    def quotient_string(self):
        return f'{self.quotient:.2f}'

    @property
    def synthese_anzeige(self):
        if not self.ergebnis_ist_vorhanden():
            return ''
        if self.synthese_pct == 0:
            return '0%'
        return f'{self.synthese_pct*100:.1f}%'

    def ergebnis_ist_vorhanden(self) -> bool:
        return self.value_csf is not None

    def berechne_quotient(self):
        if self.value_csf and self.value_serum:
            self.quotient = float(self.value_csf) / float(self.value_serum)


@dataclass(frozen=True)
class PatientData:
    """
    Class holding all necessary informations for creating a Reiber diagram

    Attributes
    ----------
    birth_date_iso : str
                     String like 2002-09-16

    albumin_serum : float
                    concentration of Albumin in serum sample
    albumin_csf : float
                  concentration of Albumin in CSF sample. Must have same unit as albumin_serum
    iga_serum : float
                concentration of IgA in serum sample
    iga_csf : float
              concentration of IgA in CSF sample. Must have same unit as iga_serum
    igg_serum : float
                concentration of IgG in serum sample
    igg_csf : float
              concentration of IgG in CSF sample. Must have same unit as igg_serum
    igm_serum : float
                concentration of IgG in serum sample
    igm_csf : float
              concentration of IgM in CSF sample. Must have same unit as igm_serum
    """
    birth_date_iso: str
    albumin_serum: float
    albumin_csf: float
    igg_serum: Optional[float] = None
    igg_csf: Optional[float] = None
    iga_serum: Optional[float] = None
    iga_csf: Optional[float] = None
    igm_serum: Optional[float] = None
    igm_csf: Optional[float] = None

    @property
    def albumin_quotient(self) -> float:
        return self.albumin_csf / self.albumin_serum

    @property
    def iga_quotient(self) -> float:
        return self.iga_csf / self.iga_serum

    @property
    def igg_quotient(self) -> float:
        return self.igg_csf / self.igg_serum

    @property
    def igm_quotient(self) -> float:
        return self.igm_csf / self.igm_serum

    def get_ig_quotient(self, name: Immunglobulin) -> float:
        if name == Immunglobulin.IGA:
            return self.iga_quotient
        elif name == Immunglobulin.IGG:
            return self.igg_quotient
        elif name == Immunglobulin.IGM:
            return self.igm_quotient
        return None

    @staticmethod
    def _only_one_value_exists(val1: Optional[float], val2: Optional[float]) -> bool:
        if val1 and not val2:
            return True
        if not val1 and val2:
            return True
        return False

    def is_valid(self) -> bool:
        if self._only_one_value_exists(self.igg_serum, self.igg_csf):
            return False
        elif self._only_one_value_exists(self.iga_serum, self.iga_csf):
            return False
        elif self._only_one_value_exists(self.igm_serum, self.igm_csf):
            return False
        return True

    @property
    def birth_date(self) -> datetime.datetime:
        return datetime.datetime.strptime(self.birth_date_iso, "%Y-%m-%d")

    @property
    def age(self) -> int:
        now = datetime.datetime.now()
        delta = now - self.birth_date
        return int(delta.days / 365.25)

    def albumin_result(self) -> AnalytResult:
        return AnalytResult('Albumin', value_serum=self.albumin_serum, value_csf=self.albumin_csf)

    def to_result_list(self) -> list[AnalytResult]:
        result_list = []
        if self.iga_serum:
            result_list.append(AnalytResult(Immunglobulin.IGA, value_serum=self.iga_serum, value_csf=self.iga_csf))
        if self.igg_serum:
            result_list.append(AnalytResult(Immunglobulin.IGG, value_serum=self.igg_serum, value_csf=self.igg_csf))
        if self.igm_serum:
            result_list.append(AnalytResult(Immunglobulin.IGM, value_serum=self.igm_serum, value_csf=self.igm_csf))
        return result_list

    @classmethod
    def from_csv_file(cls, line: str) -> 'PatientData':
        data = line.split(";")
        iga_serum = float(data[3]) if len(data) > 3 else None
        iga_csf = float(data[4]) if len(data) > 4 else None
        igg_serum = float(data[5]) if len(data) > 5 else None
        igg_csf = float(data[6]) if len(data) > 6 else None
        igm_serum = float(data[7]) if len(data) > 7 else None
        igm_csf = float(data[8]) if len(data) > 8 else None
        return cls(birth_date_iso=data[0], albumin_serum=float(data[1]), albumin_csf=float(data[2]),
                   iga_serum=iga_serum, iga_csf=iga_csf, igg_serum=igg_serum, igg_csf=igg_csf, igm_serum=igm_serum,
                   igm_csf=igm_csf)


@dataclass(frozen=True)
class AnalytIdentifikation:
    verf_liquor: str
    verf_serum: str
    code_liquor: str
    code_serum: str






def _create_orders_from_csv(file_path: str):
    with open(file_path, 'r') as csv_datei:
        lines = csv_datei.readlines()
        for line in lines:
            auftrag_list.append(PatientData.from_csv_file(line))

    #for auftrag_nr, auftrag in auftrag_map.items():
    #    auftrag.berechne_delpq()


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


def create_images(data: PatientData, out_file: str, image_type: ImageType):
    """
    Parameters
    ----------
    data : PatientData
           object containing all necessary information like age and measurement data
    out_file : str
               either an absolute path if only IgG is provided or a template like
               <dir_path>/<some_name> => <dir_path>/<some_name>_<igx>.<image_type>
    image_type : ImageType
                 file extension (see ImageType class)
    """
    # parameters of the hyperbolic curves:
    immun_globulin_param_map = dict(
        IgG=CurveParameters(Immunglobulin.IGG, 0.33, 2, 0.3, 0.93, 6, 1.7),
        IgA=CurveParameters(Immunglobulin.IGA, 0.17, 74, 1.3, 0.77, 23, 3.1),
        IgM=CurveParameters(Immunglobulin.IGM, 0.04, 442, 0.82, 0.67, 120, 7.1)
    )

    # determine age-dependent vertical line:
    q_alb = (4 + data.age/15) * 0.001

    y_max = 0.15

    # Schleife ueber alle Immunglobuline:
    results = data.to_result_list()
    for analyt_ergebnis in results:

        ig_select = analyt_ergebnis.name

        if not analyt_ergebnis.ergebnis_ist_vorhanden():
            continue

        igs_curves = ImmunglobulineCurves(immun_globulin_param_map[ig_select])

        #igs_curves.setze_synthese_pct_fuer_ergebnis(analyt_ergebnis, data.albumin_result())
        #print("synthese", analyt_ergebnis.synthese_pct)

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
            print(data.albumin_quotient, data.get_ig_quotient(ig_select))
            plt.scatter(data.albumin_quotient, data.get_ig_quotient(ig_select), c='black', marker='D', s=100)

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

        file_name, file_extension = os.path.splitext(out_file)
        if file_extension:
            file_extension = file_extension.replace(".", "")
        else:
            file_extension = image_type
            out_file = f"{out_file}_{ig_select}.{file_extension}"

        plt.savefig(out_file, format=file_extension,  bbox_inches='tight')


def create_images_for_file(file_path: str, out_file: str, image_type: ImageType):
    """
    Extracts data from a csv file with order of columns (csv file has to use ; as delimiter)
    birthdate_iso | albumin_serum | albumin_csf | iga_serum | iga_csf | igg_serum | igg_csf | igm_serum | igm_csf
    Example (just IgG provided):
    2002-09-18;49.2;1.3;;;1090.1;22.3
    Parameters
    ----------
    file_path : str
               absolute path to the csv file
    out_file : str
               either an absolute path if only IgG is provided or a template like
               <dir_path>/<some_name> => <dir_path>/<some_name>_<igx>.<image_type>
    image_type : ImageType
                 file extension (see ImageType class)
    """
    _create_orders_from_csv(file_path)
    for auftrag in auftrag_list:
        create_images(auftrag, out_file, image_type)

