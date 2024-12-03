import requests
import re
from requests.adapters import HTTPAdapter, Retry
from RAPDOR.datastructures import RAPDORData
from scipy.optimize import curve_fit
from scipy.stats import norm
import plotly.graph_objs as go
import plotly.express as px
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy import integrate

upper_bound = 0
lower_bound = 20


def estimate_initial_parameters(num_components):
    # Estimate amplitudes
    guess = []
    means = np.linspace(0, 20, num_components)
    for i in range(1, num_components + 1):
        guess.append(1 / num_components)
        guess.append(means[i - 1])
        guess.append(1)

    return guess


def multi_gaussian(x, *params):
    num_gaussians = len(params) // 3
    result = np.zeros_like(x)
    for i in range(num_gaussians):
        amp, mean, std = params[i * 3:i * 3 + 3]
        result += amp * norm.pdf(x, loc=mean, scale=std)
    total_amplitude = np.sum(params[::3])  # Sum all amplitudes
    result /= total_amplitude
    return result

def genStaticFct(params):
    def mg(x):
        num_gaussians = len(params) // 3
        result = np.zeros_like(x)
        for i in range(num_gaussians):
            amp, mean, std = params[i * 3:i * 3 + 3]
            result += amp * norm.pdf(x, loc=mean, scale=std)
        total_amplitude = np.sum(params[::3])  # Sum all amplitudes
        result /= total_amplitude
        return result
    return mg

def calculate_r_squared(y_true, y_pred):
    """Calculate R-squared given the observed and predicted values."""
    rss = np.sum((y_true - y_pred)**2)  # Residual sum of squares
    tss = np.sum((y_true - np.mean(y_true))**2)  # Total sum of squares
    r_squared = 1 - (rss / tss)
    return r_squared


def multicore_fit_wrapper(x, y):
    best_bic = np.inf
    best_params = None
    best = None
    num_p = None
    for p in range(1, 6):

        initial_guess = estimate_initial_parameters(p)
        bounds = ([0, 0, 0] * p, [1, 20, np.inf] * p)
        try:

            params, cov = curve_fit(multi_gaussian, x, y, p0=initial_guess, bounds=bounds)
            fit = multi_gaussian(x, *params)
            residuals = y - fit
            #integral, _ = integrate.quad(multi_gaussian, upper_bound, lower_bound, args=params, points=[0,20] )
            #print(f"The integral of the function is approximately: {integral}")
            sse = np.sum(residuals ** 2)
            n = len(y)
            k = len(params)
            bic = n * np.log(sse / n) + k * np.log(n)
            if bic < best_bic:
                best_bic = bic
                best_params = params
                best = fit
                num_p = p
        except (RuntimeError, ValueError):
            pass
    r_squared = calculate_r_squared(y, best) if best is not None else None

    return r_squared, num_p, best_params


def determineFits(file, nr_cores=7):
    with open(file) as handle:
        data = RAPDORData.from_json(handle.read())
    i = data.state.kernel_size // 2
    x = data.fractions[i:-i].astype(float)

    sample = []
    treatments = []
    calls = []
    proteins=[]
    for idx, protein in enumerate(data.norm_array):
        for tidx, indices in enumerate(data.indices):
            treatment = data.treatment_levels[tidx]
            for inner_idx in indices:
                y = data.norm_array[idx][inner_idx]
                protein = data.df.loc[idx]["old_locus_tag"]
                calls.append((x, y))
                treatments.append(treatment)
                proteins.append(protein)
                sample.append(data.internal_design_matrix.iloc[inner_idx]["Replicate"])
    with Pool(processes=nr_cores) as pool:
        result = pool.starmap(multicore_fit_wrapper, calls)
    r_squareds, ps, best_fits = zip(*result)
    df = pd.DataFrame(
        {
            "old_locus_tag": proteins,
            "Treatment": treatments,
            "Replicate": sample,
            "fit": best_fits,
            "R-squared": r_squareds,
            "num_gauss": ps
        }
    )
    df = df.sort_values(by="R-squared", ascending=True)
    return df


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
    url = f"https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Cgene_oln%2Corganism_name%2Clength%2Cgo_id%2Cft_transmem%2Cft_intramem%2Cxref_string&format=tsv&query=%28%28taxonomy_id%3A{tax_id}%29%29&size=500"
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


def enrichment_plot_from_cp_table(df, mode="scatter",  colors=("Red", "Blue")):
    if len(df) == 0:
        return empty_figure()
    def df_div(l):
        return int(l[0]) / int(l[1])
    df["GeneRatio"] = df["GeneRatio"].str.split("/").apply(df_div)
    if mode == "scatter":
        fig = px.scatter(
            df,
            x="GeneRatio",
            y="Description",
            symbol="ONTOLOGY",
            color="p.adjust",
            template="plotly_white",

        )
        fig.update_traces(marker=dict(size=15))
    elif mode == "bar":
        fig = px.bar(
            df,
            x="GeneRatio",
            y="Description",
            color="p.adjust",
            template="plotly_white",
            color_continuous_scale=[colors[0], colors[1]]

        )
    else:
        raise ValueError(f"mode: {mode} is not valid")


    fig.update_layout(
        coloraxis_colorbar=dict(
            yanchor="top",
            y=0.7,
            len=0.7,
            x=1,
            ticks="outside"
        ),
        legend=dict(x=1),
        yaxis=dict(tickmode="linear", type="category", dtick=1)
    )
    return fig

if __name__ == '__main__':

    #download_organism_go_terms(1111708, outfile="../Data/GOAnno.tsv")
    download_organism_go_terms(9606, outfile="../Data/GOAnnoHuman.tsv")