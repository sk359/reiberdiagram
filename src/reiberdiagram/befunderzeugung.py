from jinja2 import Environment, FileSystemLoader

env = Environment(loader=FileSystemLoader(self.settings.project_root))
template = env.get_template(f"src/templates/weiblich_einzel.html")

html_string = template.render(date=datum_heute, auftrag_nr=befund_data.auftrag_nr,
                              befund_datum=befund_data.befund_datum,
                              freigabedatum=befund_data.freigabe_datum,
                              patientin=befund_data.patientin, einsender=befund_data.einsender,
                              unterschriften_ordner=untr,
                              noch_nicht_freigegeben=befund_data.ist_nicht_freigegeben())
# html_string = html_string.encode('utf8').decode('utf8')
temp_file = "temp.html"
with open(temp_file, "w") as tfh:
    tfh.write(html_string)

htmldoc = HTML(filename=temp_file, base_url="")
pdf = htmldoc.write_pdf(target=self.pdf_file_path)
