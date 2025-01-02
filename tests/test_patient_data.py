import pytest
import datetime

from reiberdiagram import PatientData, Immunglobulin, InvalidData

class TestPatientDataSuite:

    def test_class(self):
        now = datetime.datetime.now()
        born_days_ago = int(365*29.5)
        birth_date = now - datetime.timedelta(days=born_days_ago)
        pat_data = PatientData(birth_date_iso=birth_date.strftime("%Y-%m-%d"), albumin_serum=100, albumin_csf=0.5,
                               igg_serum=200, igg_csf=2.0,
                               iga_serum=50, iga_csf=0.1,
                               igm_serum=1000, igm_csf=0.1)
        assert pat_data.albumin_quotient == 0.005
        assert pat_data.igg_quotient == 0.01
        assert pat_data.iga_quotient == 0.002
        assert pat_data.igm_quotient == 0.0001
        assert pat_data.get_ig_quotient(Immunglobulin.IGA) == 0.002
        assert pat_data.age == 29

    def test_invalid_data(self):
        pat_data = PatientData(birth_date_iso="2002-01-01", albumin_serum=100, albumin_csf=0.5,
                               igg_serum=200)
        with pytest.raises(InvalidData):
            pat_data.raise_if_invalid()

    def test_to_list(self):
        pat_data = PatientData(birth_date_iso="2002-01-01", albumin_serum=100, albumin_csf=0.5,
                               igg_serum=200, igg_csf=2.0,
                               iga_serum=50, iga_csf=0.1,
                               igm_serum=1000, igm_csf=0.1)

        result_list = pat_data.to_result_list()
        assert len(result_list) == 3
        assert result_list[0].name == Immunglobulin.IGA