import os
import pytest

from reiberdiagram import create_images, PatientData, ImageType

test_dir = os.path.dirname(os.path.realpath(__file__))
test_file_path = f"{test_dir}/test"
test_file = f"{test_file_path}_IgG.png"
base_file = f"{test_dir}/baseline_IgG.png"


@pytest.fixture(autouse=True)
def run_around_tests():
    yield
    try:
        os.remove(test_file)
    except OSError:
        pass


def test_image_creation():
    pat = PatientData(birth_date_iso="2022-09-13", albumin_serum=1000, albumin_csf=10, igg_serum=200, igg_csf=1)
    create_images(data=pat, out_file=test_file_path, image_type=ImageType.PNG)
    assert os.path.exists(test_file)
    assert os.path.getsize(test_file) == os.path.getsize(base_file)