![Codecov](https://img.shields.io/codecov/c/github/sk359/reiberdiagram)


# Reiberdiagram

Python module to create an image file for each measured immunoglobulin.

The Reiberdiagram was invented by Hansotto Reiber to visualize reasons for unusual high concentrations of proteins
in the Cerebrospinal fluid (CSF).[^1][^2] This may be due to problems with the blood brain barrier or antibody
synthesis in the CSF. The diagram consist of two main curves enclosing the normal range, some additional curves
labelled 20%, 40% etc., an age-dependent vertical separation line and a mark that highlights the position of the 
current measurement. The x-axis is the albumin quotient and the y-axis the immunoglobulin quotient.

Images can be created for IgA, IgG and IgM

### Usage:

You can either evaluate a csv file (see example in test folder) or a single patient. Example:

```
test_csv_file:
id;dob;alb_serum;alb_csf;iga_serum;iga_csf;igg_serum;igg_csf;igm_serum;igm_csf
001;2002-09-18;49.2;1.3;;;1090.1;22.3
002;1976-11-19;49.2;1.3;100.2;13.4;;;345.6;33.1
```

```python
from reiberdiagram import create_images_for_file, ImageType
create_images_for_file(test_csv_file, out_file=test_file_path, image_type=ImageType.PNG)
```

```python
from reiberschema import create_images, PatientData, ImageType
pat = PatientData(csv_row_id="001", birth_date_iso="2022-09-13", albumin_serum=1000, albumin_csf=10, igg_serum=133.5, igg_csf=1.5)
create_images(data=pat, out_file="test", image_type=ImageType.PNG)
```

### Example image:

![Example diagram for IgG](/tests/baseline_IgG.png)

[^1]: Reiber H (1994). The hyperbolic function: a mathematical solution of the protein flux/CSF flow model 
for blood-CSF barrier function J Neurol Sci 126:243-245.

[^2]: Reiber H (1994). Flow rate of cerebrospinal fluid (CSF)- a concept common to normal blood-CSF barrier function 
and to dysfunction in neurological diseases. J Neurol Sci 122:189-203
