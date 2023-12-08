import requests
import re
from requests.adapters import HTTPAdapter, Retry



def download_organism_go_terms(tax_id, outfile):
    re_next_link = re.compile(r'<(.+)>; rel="next"')

    def get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = get_next_link(response.headers)
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    #url = f"https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Cgene_oln%2Corganism_name%2Clength%2Cgo_id%2Cgo&format=tsv&query=%28%28taxonomy_id%3A{tax_id}08%29%29&size=500"
    url = f"https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Cgene_oln%2Corganism_name%2Clength%2Cgo_id%2Cft_transmem%2Cft_intramem&format=tsv&query=%28%28taxonomy_id%3A{tax_id}08%29%29&size=500"
    progress = 0
    print("Starting to Download GO Terms")
    with open(outfile, "w") as f:
        for batch, total in get_batch(url):
            lines = batch.text.splitlines()
            if not progress:
                print(lines[0], file=f)
            for line in lines[1:]:
                print(line, file=f)
            progress += len(lines[1:])
            print(f'Downloaded {progress} / {total}')
    return 0

if __name__ == '__main__':

    download_organism_go_terms(11117, outfile="../Data/GOAnno.tsv")