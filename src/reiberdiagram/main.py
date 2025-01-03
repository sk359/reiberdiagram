import os
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional
import math
import datetime
from dataclasses import dataclass
import bisect

from constants import ImageType, Immunglobulin

base_path_reiberschema = os.path.abspath(os.path.join(os.path.dirname(__file__)))

order_list = []


class InvalidData(Exception):
    pass


@dataclass(frozen=True)
class PatientData:
    """
    Class holding all necessary information for creating a Reiber diagram

    Attributes
    ----------
    birth_date_iso : str
                     String like 2002-09-16
    csv_row_id : str
                 will be part of file name, can be any string (no spaces)
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
    csv_row_id: str
    albumin_serum: float
    albumin_csf: float
    igg_serum: Optional[float] = None
    igg_csf: Optional[float] = None
    iga_serum: Optional[float] = None
    iga_csf: Optional[float] = None
    igm_serum: Optional[float] = None
    igm_csf: Optional[float] = None

    def __repr__(self) -> str:
        return str(vars(self))

    def has_immunoglobulin_result(self, name: Immunglobulin) -> bool:
        if name == Immunglobulin.IGA:
            return self.iga_serum is not None
        elif name == Immunglobulin.IGG:
            return self.igg_serum is not None
        elif name == Immunglobulin.IGM:
            return self.igm_serum is not None
        return False

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

    def raise_if_invalid(self):
        if not self.is_valid():
            raise InvalidData("Invalid data provided (probably missing serum or csf measurement)")

    @property
    def birth_date(self) -> datetime.datetime:
        return datetime.datetime.strptime(self.birth_date_iso, "%Y-%m-%d")

    @property
    def age(self) -> int:
        now = datetime.datetime.now()
        delta = now - self.birth_date
        return int(delta.days / 365.25)

    @staticmethod
    def string_to_float_or_none(value_list: list[str], index: int) -> Optional[float]:
        try:
            float_number = float(value_list[index])
            return float_number
        except (ValueError, IndexError):
            return None

    @classmethod
    def from_csv_file(cls, line: str) -> 'PatientData':
        data = line.split(";")
        iga_serum = PatientData.string_to_float_or_none(data, 4)
        iga_csf = PatientData.string_to_float_or_none(data, 5)
        igg_serum = PatientData.string_to_float_or_none(data, 6)
        igg_csf = PatientData.string_to_float_or_none(data, 7)
        igm_serum = PatientData.string_to_float_or_none(data, 8)
        igm_csf = PatientData.string_to_float_or_none(data, 9)
        return cls(csv_row_id=data[0], birth_date_iso=data[1], albumin_serum=float(data[2]), albumin_csf=float(data[3]),
                   iga_serum=iga_serum, iga_csf=iga_csf, igg_serum=igg_serum, igg_csf=igg_csf, igm_serum=igm_serum,
                   igm_csf=igm_csf)


@dataclass(frozen=True)
class AnalytIdentifikation:
    verf_liquor: str
    verf_serum: str
    code_liquor: str
    code_serum: str


class DiagramDimension:
    x_values = [x * 0.0001 for x in range(10, 1000)]
    y_max = 0.15


class CurveParameters:

    """
    Parameters for a given Immunoglobulin based on
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
    Holds list of x and y values that are used inside the diagram for the curves
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
               <dir_path>/<some_name> => <dir_path>/<some_name>_<row_id>_<igx>.<image_type>
    image_type : ImageType
                 file extension (see ImageType class)
    """
    data.raise_if_invalid()
    out_file_new = out_file

    # parameters of the hyperbolic curves:
    immun_globulin_param_map = dict(
        IgG=CurveParameters(Immunglobulin.IGG, 0.33, 2, 0.3, 0.93, 6, 1.7),
        IgA=CurveParameters(Immunglobulin.IGA, 0.17, 74, 1.3, 0.77, 23, 3.1),
        IgM=CurveParameters(Immunglobulin.IGM, 0.04, 442, 0.82, 0.67, 120, 7.1)
    )

    # determine age-dependent vertical line:
    q_alb = (4 + data.age/15) * 0.001

    y_max = 0.15

    ig_list = [Immunglobulin.IGA, Immunglobulin.IGG, Immunglobulin.IGM]
    for immunglobulin in ig_list:

        if not data.has_immunoglobulin_result(immunglobulin):
            continue

        igs_curves = ImmunglobulineCurves(immun_globulin_param_map[immunglobulin])

        fig = plt.figure()
        plt.rc('axes', axisbelow=True)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(igs_curves.x, igs_curves.lower, color=(0, 0, 0))
        ax.set_xscale('log')
        ax.set_yscale('log')
        # Reference section (upper and lower thick curve):
        ax.plot(igs_curves.x, igs_curves.lower, color=(0, 0, 0), linewidth=3)
        ax.plot(igs_curves.x, igs_curves.upper, color=(0, 0, 1), linewidth=3)

        # Percentage curves intrathekale fraction:
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

        # Highlight measurement point:
        plt.scatter(data.albumin_quotient, data.get_ig_quotient(immunglobulin), c='black', marker='D', s=100)

        # Axis:
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
            out_file_new = f"{out_file}_{data.csv_row_id}_{immunglobulin}.{file_extension}"

        plt.savefig(out_file_new, format=file_extension,  bbox_inches='tight')


def _create_orders_from_csv(file_path: str):
    def is_header_line(csv_line: str) -> bool:
        # returns True if the seconds column does not hold an isoformat date string
        cols = line.split(";")
        try:
            _ = datetime.datetime.strptime(cols[1], "%Y-%m-%d")
            return False
        except Exception:
            return True

    with open(file_path, 'r') as csv_datei:
        lines = csv_datei.readlines()
        for row_count, line in enumerate(lines):
            if row_count == 0 and is_header_line(line):
                continue
            order_list.append(PatientData.from_csv_file(line))


def create_images_for_file(file_path: str, out_file: str, image_type: ImageType):
    """
    Extracts data from a csv file with order of columns (csv file has to use ; as delimiter)
    id | birthdate_iso | albumin_serum | albumin_csf | iga_serum | iga_csf | igg_serum | igg_csf | igm_serum | igm_csf
    Example (just IgG provided):
    001;2002-09-18;49.2;1.3;;;1090.1;22.3

    The ID can have any format. A header row can be provided in the csv file but is optional

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
    for order in order_list:
        create_images(order, out_file, image_type)
