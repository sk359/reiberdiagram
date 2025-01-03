import os
import pytest

from reiberdiagram import create_images, create_images_for_file, PatientData, ImageType

test_dir = os.path.dirname(os.path.realpath(__file__))
test_file_path = f"{test_dir}/test"
test_file = f"{test_file_path}_01_IgG.png"
base_file = f"{test_dir}/baseline_IgG.png"
test_csv_file = f"{test_dir}/testfile.csv"


@pytest.fixture(autouse=True)
def remove_test_images():
    yield
    try:
        for file in os.listdir(test_dir):
            if file.endswith(".png"):
                file_abs_path = os.path.join(test_dir, file)
                if file_abs_path == base_file:
                    continue
                os.remove(file_abs_path)
    except OSError:
        pass

def get_all_files_by_ending(ending: str):
    file_list = []
    for file in os.listdir(test_dir):
        if file.endswith(ending):
            file_list.append(file)
    return file_list


def test_image_creation():
    pat = PatientData(csv_row_id="01", birth_date_iso="2022-09-13", albumin_serum=1000, albumin_csf=10, igg_serum=200,
                      igg_csf=1)
    create_images(data=pat, out_file=test_file_path, image_type=ImageType.PNG)
    assert os.path.exists(test_file)
    assert os.path.getsize(test_file) == os.path.getsize(base_file)

def test_image_creation_from_csv_file():
    create_images_for_file(test_csv_file, out_file=test_file_path, image_type=ImageType.PNG)
    file_list = get_all_files_by_ending(".png")
    expected_numer_of_images = 4 # 3 + test baseline image
    assert len(file_list) == expected_numer_of_images
    assert "test_001_IgG.png" in file_list
    assert "test_002_IgM.png" in file_list