import pandas as pd
from jinja2 import PackageLoader, Environment, FileSystemLoader

def generate_html(knead_result, kraken_result,prokka_result, megahit_result):
    env = Environment(loader=FileSystemLoader('/cygene/work/00.test/pipeline/meta_genomics_pipeline/database/templates'))
    template = env.get_template('index.html')
    with open("./summary_report/result.html","w+") as fout:
        html_content = template.render(knead_result=knead_result, 
                kraken_result=kraken_result,
                prokka_result=prokka_result, 
                megahit_result = megahit_result)
        fout.write(html_content)

def deal_data(data):
    data = pd.read_table(data,sep="\t")
    # data = data.sort_values(by='Sample')
    data_result = data.to_html(index=False)
    return data_result

def copy_dependency():
    sys.(cp  summary_dir)
    # copy dependency to html result file.
    if :
        sys.system("cp templets_folder summary_dir".format())
    s

knead_result = deal_data("Summary.kneadata.results.xls")
kraken_result = deal_data("Summary.kraken2.results.xls")
#humann_result = deal_data("Summary.humann2.results.xls")
prokka_result = deal_data("Summary.prokka.results.xls")
megahit_result = deal_data("Summary.megahit.quest.results.xls")


generate_html(knead_result, kraken_result, prokka_result, megahit_result)
